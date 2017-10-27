# pylint: disable=too-many-arguments,invalid-name,too-many-locals,too-many-branches
"""
Defines functions and classes to start the processing of Vesuvio data.
The main entry point that most users should care about is fit_tof().
"""
from __future__ import (absolute_import, division, print_function)
from six import iteritems

import copy
import re
import numpy as np

from mantid import mtd
from mantid.api import (AlgorithmManager, AnalysisDataService, MatrixWorkspace, WorkspaceFactory, TextAxis)
from mantid.kernel import MaterialBuilder
from vesuvio.instrument import VESUVIO

from mantid.simpleapi (CropWorkspace, ConjoinWorkspaces, DeleteWorkspace, GroupWorkspaces, LoadVesuvio, Rebin,
                       RenameWorkspace, UnGroupWorkspace, VesuvioCorrections, VesuvioTOFFit)


# --------------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------------

def fit_tof(runs, flags, iterations=1, convergence_threshold=None):
    vesuvio_loader = VesuvioLoadHelper(flags['diff_mode'], flags['fit_mode'],
                                       flags['ip_file'], flags['bin_parameters'])
    vesuvio_input = VesuvioTOFFitInput(runs, flags['container_runs'],
                                       flags['spectra'], vesuvio_loader)
    ms_helper = None
    if flags.get('ms_enabled', True):
        hydrogen_constraints = flags['ms_flags'].pop("HydrogenConstraints", {})
        ms_helper = VesuvioMSHelper(**flags['ms_flags'])
        ms_helper.add_hydrogen_constraints(hydrogen_constraints)

    fit_helper = VesuvioTOFFitHelper(_create_background_str(flags.get('background', None)),
                                     _create_intensity_constraint_str(flags['intensity_constraints']),
                                     _create_user_defined_ties_str(flags['masses']),
                                     flags.get('max_fit_iterations', 5000),
                                     flags['fit_minimizer'])

    corrections_helper = VesuvioCorrectionsHelper(flags.get('gamma_correct', False),
                                                  flags.get('ms_enabled', True),
                                                  flags.get('fixed_gamma_scaling', 0.0),
                                                  flags.get('fixed_container_scaling', 0.0),
                                                  flags['intensity_constraints'])

    mass_profile_collection = MassProfileCollection(flags['masses'])
    vesuvio_fit_routine = VesuvioTOFFitRoutine(ms_helper, fit_helper, corrections_helper,
                                               mass_profile_collection, flags['fit_mode'])
    vesuvio_output, exit_iteration  = vesuvio_fit_routine(vesuvio_input, iterations, convergence_threshold,
                                                          flags.get('output_verbose_corrections', False))
    return (vesuvio_output.result_workspace(), vesuvio_output.prefit_parameters_workspace(),
            vesuvio_output.parameters_workspace(), exit_iteration)


# -----------------------------------------------------------------------------------------

class VesuvioTOFFitInput(object):
    """
    A helper class for loading and storing the specified spectra of the specified
    input sample and container runs, using the given loader.

    Attributes:
        sample_runs         The sample runs to load and store.
        container_runs      The container runs to load and store.
        sample_data         The loaded sample runs as a workspace.
        container_data      The loaded container runs as a workspace
        spectra             The spectra to load.
        _back_scattering    True if back scattering spectra are being
                            loaded, false otherwise.
    """

    def __init__(self, sample_runs, container_runs, spectra, loader):
        # Load sample and container runs
        self.sample_runs, self.sample_data = \
            self._try_load_data(sample_runs, spectra, loader,
                                "Unable to load Sample Runs:")
        mtd.addOrReplace(self._tof_workspace_suffix(sample_runs, spectra),
                         self.sample_data)

        if container_runs is not None:
            self.container_runs, self.container_data = \
                self._try_load_data(container_runs, spectra, loader,
                                    "Unable to load Container Runs:")
            mtd.addOrReplace(self._tof_workspace_suffix(container_runs, spectra),
                             self.container_data)
        else:
            self.container_runs, self.container_data = None, None

        self.spectra = spectra
        self._back_scattering = self._is_back_scattering_spectra(spectra)

    def using_back_scattering_spectra(self):
        """
        :return: True if back-scattering spectra is being used as input,
                 False otherwise.
        """
        return self._back_scattering

    def _is_back_scattering_spectra(self, spectra):
        """
        Checks whether the specified spectra denotes back scattering spectra.

        :param spectra: The spectra to check.
        :return:        True if the specified spectra denotes back scattering
                        spectra.
        """
        if isinstance(spectra, str):
            return spectra == 'backward'
        else:
            try:
                first_spec = int(spectra.split("-")[0])
                back_banks = VESUVIO().backward_banks
                return any([lower <= first_spec <= upper for lower, upper in back_banks])
            except:
                raise RuntimeError("Invalid value given for spectrum range: Range must "
                                   "either be 'forward', 'backward' or specified with "
                                   "the syntax 'a-b'.")

    def _try_load_data(self, runs, spectra, loader, error_msg):

        try:
            return self._load_data(runs, spectra, loader)
        except Exception as exc:
            raise RuntimeError(error_msg + "\n" + str(exc))

    def _load_data(self, runs, spectra, loader):
        """
        Loads the specified spectra, of the specified runs, using the
        specified loader.

        :param runs:    The runs to load - either an already loaded
                        runs workspace or an input string of runs.
        :param spectra: The spectra of the runs to load.
        :param loader:  The loader to use.
        :return:        A tuple of the name of the input workspace
                        or input string and the loaded output workspace.
        """
        if isinstance(runs, MatrixWorkspace):
            return runs.getName(), runs
        else:
            return runs, loader(runs, spectra)

    def _tof_workspace_suffix(self, runs, spectra):
        """
        Creates and returns the suffix for the loaded TOF workspaces.

        :param runs:    The runs being loaded.
        :param spectra: The spectra being loaded.
        :return:        A suffix for the TOF workspaces being loaded.
        """
        return runs + "_" + spectra + "_tof"


# -----------------------------------------------------------------------------------------

class VesuvioTOFFitRoutine(object):
    """
    A class for executing the the Vesuvio TOF Fit Routine from a Vesuvio Driver Script.

    Attributes:
        _ms_helper                   A helper object for multiple scattering parameters.
        _fit_helper                  A helper object for computing a VesuvioTOFFit.
        _corrections_helper          A helper object for computing VesuvioCorrections.
        _mass_profile_collection     An object for storing and manipulating mass values
                                     and profiles.
        _fit_mode                    The fit mode to use in the fitting routine.
    """

    def __init__(self, ms_helper, fit_helper, corrections_helper, mass_profile_collection, fit_mode):
        self._ms_helper = ms_helper
        self._fit_helper = fit_helper
        self._corrections_helper = corrections_helper
        self._mass_profile_collection = mass_profile_collection
        self._fit_mode = fit_mode

    def __call__(self, vesuvio_input, iterations, convergence_threshold, verbose_output=False):
        if iterations < 1:
            raise ValueError('Must perform at least one iteration')
        # Creation of a fit routine iteration
        tof_iteration = VesuvioTOFFitRoutineIteration(self._ms_helper, self._fit_helper,
                                                      self._corrections_helper,
                                                      self._mass_profile_collection,
                                                      self._fit_mode)

        back_scattering = vesuvio_input.using_back_scattering_spectra()
        update_filter = ignore_hydrogen_filter if back_scattering else None
        exit_iteration = 0
        vesuvio_output = VesuvioTOFFitOutput(vesuvio_input, verbose_output)

        for iteration in range(1, iterations + 1):
            # Update the mass profiles using the previous result if it exists
            if vesuvio_output.contains_fit() is not None:
                self._mass_profile_collection.update_profiles_from_workspace(previous_results[2],
                                                                             update_filter)
            print("=== Iteration {0} out of a possible {1}".format(iteration, iterations))
            vesuvio_output = tof_iteration(vesuvio_input, vesuvio_output, iteration, verbose_output)
            exit_iteration += 1
            # Check whether the change in the cost function between the result and
            # previous results is smaller than the convergence threshold.
            if vesuvio_output.contains_fit() and convergence_threshold is not None:
                cost_function_change = vesuvio_output.change_in_cost_function()
                print("Cost function change: {0}".format(cost_function_change))

                if cost_function_change <= convergence_threshold:
                    print("Stopped at iteration {0} due to minimal change in cost function"
                          .format(exit_iteration))
                    return vesuvio_output, exit_iteration
            vesuvio_output.increment_iteration()
        return vesuvio_output, exit_iteration


# ------------------------------------------------------------------------------------------------------

class VesuvioTOFFitRoutineIteration(object):
    """
    A class for executing a single iteration of the Vesuvio TOF Fit Routine, from a
    Vesuvio Driver Script.

    Attributes:
        _ms_helper                   A helper object for multiple scattering parameters.
        _fit_helper                  A helper object for computing a VesuvioTOFFit.
        _corrections_helper          A helper object for computing VesuvioCorrections.
        _mass_profile_collection     An object for storing and manipulating mass values
                                     and profiles.
        _fit_mode                    The fit mode to use in the fitting routine.
    """

    def __init__(self, ms_helper, fit_helper, corrections_helper, mass_profile_collection, fit_mode):
        self._ms_helper = ms_helper
        self._fit_helper = fit_helper
        self._corrections_helper = corrections_helper
        self._mass_profile_collection = mass_profile_collection

    def __call__(self, vesuvio_input, vesuvio_output, iteration, verbose_output=False):
        # Retrieve all input in local variables
        sample_runs = str(vesuvio_input.sample_runs)
        sample_data = vesuvio_input.sample_data
        container_data = vesuvio_input.container_data
        num_spectra = sample_data.getNumberHistograms()

        fit_filter = ignore_hydrogen_filter if vesuvio_input.using_back_scattering_spectra() else None
        all_mass_values = self._mass_profile_collection.masses()
        fit_mass_values = self._mass_profile_collection.masses(fit_filter)

        for index in range(num_spectra):
            all_profiles = self._mass_profile_collection.profile_functions(index)
            fit_profiles = self._mass_profile_collection.profile_functions(index, fit_filter)

            # Calculate pre-fit to retrieve parameter approximations for corrections
            suffix = self._create_fit_workspace_suffix(index, sample_data, self._fit_mode, iteration)
            fit_ws_name = "__vesuvio_corrections_fit"
            pre_params_ws_name = sample_runs + "_params_pre_correction" + suffix
            self._fit_helper(InputWorkspace=sample_data,
                             WorkspaceIndex=index,
                             Masses=fit_mass_values,
                             MassProfiles=fit_profiles,
                             OutputWorkspace=fit_ws_name,
                             FitParameters=pre_params_ws_name)
            DeleteWorkspace(fit_ws_name)

            # Calculate and apply vesuvio corrections
            linear_corrections_fit_params_name = sample_runs + "_correction_fit_scale" + suffix
            corrected_data_name = sample_runs + "_tof_corrected" + suffix
            corrections_args = {'FitParameters': pre_params_ws_name}

            if self._ms_helper is not None:
                corrections_args.update(self._ms_helper.to_dict())
            if verbose_output:
                corrections_args['CorrectionWorkspaces'] = sample_runs + "_correction" + suffix
                corrections_args['CorrectedWorkspaces'] = sample_runs + "_corrected" + suffix
            if container_data is not None:
                corrections_args['ContainerWorkspace'] = container_data

            self._corrections_helper(InputWorkspace=sample_data,
                                     WorkspaceIndex=index,
                                     Masses=all_mass_values,
                                     MassProfiles=all_profiles,
                                     MassIndexToSymbolMap=corrected_data_name,
                                     OutputWorkspace=corrected_data_name,
                                     LinearFitResult=linear_corrections_fit_params_name,
                                     **corrections_args)
            # Calculate final fit
            fit_ws_name = sample_runs + "_data" + suffix
            params_ws_name = sample_runs + "_params" + suffix
            fit_result = self._fit_helper(InputWorkspace=corrected_data_name,
                                          WorkspaceIndex=0,
                                          Masses=fit_mass_values,
                                          MassProfiles=fit_profiles,
                                          OutputWorkspace=fit_ws_name,
                                          FitParameters=params_ws_name)
            DeleteWorkspace(corrected_data_name)

            # Update output with results from fit
            vesuvio_output.add_chi2_value(fit_result[-1])
            vesuvio_output.set_fit_workspace_name(fit_ws_name)
            vesuvio_output.update_prefit_parameters_workspace(mtd[pre_params_ws_name], index)
            vesuvio_output.update_parameters_workspace(mtd[params_ws_name], index)
            vesuvio_output.update_fit_workspace(mtd[linear_corrections_fit_params_name], index)

            # Delete parameter tables after use
            DeleteWorkspace(pre_params_ws_name)
            DeleteWorkspace(params_ws_name)
            DeleteWorkspace(linear_corrections_fit_params_name)

            # Add output to Analysis Data Service
            vesuvio_output.add_output_to_ads(index, mtd[corrections_args["CorrectionWorkspaces"]],
                                             mtd[corrections_args["CorrectedWorkspaces"]])
        return vesuvio_output

    def _create_fit_workspace_suffix(self, index, tof_data, fit_mode, spectra, iteration=None):
        if fit_mode == "bank":
            suffix = "_" + spectra + "_bank_" + str(index + 1)
        else:
            spectrum = tof_data.getSpectrum(index)
            suffix = "_spectrum_" + str(spectrum.getSpectrumNo())

        if iteration is not None:
            suffix += "_iteration_" + str(iteration)

        return suffix


# ------------------------------------------------------------------------------------------------------

class VesuvioLoadHelper(object):
    """
    A helper class for loading Vesuvio data from the input of a user script.

    Attributes:
        _diff_mode       The difference mode to load runs with.
        _fit_mode        The fit mode to load runs with.
        _param_file      The instrument parameter file for loading.
        _rebin_params    The parameters to use for rebinning loaded data.
                         If none, loaded data is not rebinned.
    """

    def __init__(self, diff_mode, fit_mode, param_file, rebin_params=None):
        self._diff_mode = self._parse_diff_mode(diff_mode)
        self._fit_mode = fit_mode
        self._param_file = param_file
        self._rebin_params = rebin_params
        self._instrument = VESUVIO()

    def __call__(self, runs, spectra):
        return self.load_and_crop_runs(runs, spectra)

    def load_and_crop_runs(self, runs, spectra):
        sum_spectra = (self._fit_mode == 'bank')
        loaded = self.load_runs(runs, self._parse_spectra(spectra), sum_spectra)
        cropped = self._crop_workspace(loaded)

        if self._rebin_params is not None:
            return self._rebin_workspace(cropped)
        else:
            return cropped

    def _parse_diff_mode(self, diff_mode):
        if diff_mode == "single":
            return "SingleDifference"
        elif diff_mode == "double":
            return "DoubleDifference"
        else:
            return diff_mode

    def _parse_spectra(self, spectra):
        if self._fit_mode == 'bank':
            return self._parse_spectra_bank(spectra)
        else:
            if spectra == "forward":
                return "{0}-{1}".format(*self._instrument.forward_spectra)
            elif spectra == "backward":
                return "{0}-{1}".format(*self._instrument.backward_spectra)
            else:
                return spectra

    def _parse_spectra_bank(self, spectra_bank):
        if spectra_bank == "forward":
            bank_ranges = self._instrument.forward_banks
        elif spectra_bank == "backward":
            bank_ranges = self._instrument.backward_banks
        else:
            raise ValueError("Fitting by bank requires selecting either 'forward' or 'backward' "
                             "for the spectra to load")
        bank_ranges = ["{0}-{1}".format(x, y) for x, y in bank_ranges]
        return ";".join(bank_ranges)

    def load_runs(self, runs, spectra, sum_spectra=False):
        return LoadVesuvio(Filename=runs,
                           Mode=self._diff_mode,
                           InstrumentParFile=self._param_file,
                           SpectrumList=spectra,
                           SumSpectra=sum_spectra,
                           OutputWorkspace="loaded",
                           StoreInADS=False)

    def _crop_workspace(self, workspace):
        return CropWorkspace(InputWorkspace=workspace,
                             XMin=self._instrument.tof_range[0],
                             XMax=self._instrument.tof_range[1],
                             OutputWorkspace="cropped",
                             StoreInADS=False)

    def _rebin_workspace(self, workspace):
        return Rebin(InputWorkspace=workspace,
                     Params=self._rebin_params,
                     OutputWorkspace="rebinned",
                     StoreInADS=False)


# ------------------------------------------------------------------------------------------------------

class VesuvioMSHelper(object):
    """
    A helper class for storing and manipulating the multiple scattering paramaters of
    the Vesuvio TOF Fit Routine.

    Attributes:
        _beam_radius        The radius of the neutron beam
        _sample_height      The height of the sample
        _sample_width       The width of the sample
        _sample_depth       The depth of the sample
        _sample_density     The density of the sample
        _seed               The seed to use in generating random
                            scattering events
        _num_scatters       The number of times a neutron will be scattered
        _num_runs           The number of runs
        _num_events         The number of scattering events to simulate
        _smooth_neighbours
    """

    def __init__(self, BeamRadius=2.5, SampleHeight=5.0, SampleWidth=5.0, SampleDepth=5.0,
                 SampleDensity=1.0, Seed=123456789, NumScatters=3, NumRuns=10, NumEvents=50000,
                 SmoothNeighbours=3):
        self._beam_radius = BeamRadius
        self._sample_height = SampleHeight
        self._sample_width = SampleWidth
        self._sample_depth = SampleDepth
        self._sample_density = SampleDensity
        self._seed = Seed
        self._num_scatters = NumScatters
        self._num_runs = NumRuns
        self._num_events = NumEvents
        self._smooth_neighbours = SmoothNeighbours
        # Setup hydrogen constraints
        def parser(constraint):
            return self._parse_hydrogen_constraint(constraint)
        self._hydrogen_constraints = VesuvioConstraints("HydrogenConstraints", parser)

    def add_hydrogen_constraints(self, constraints):
        self._hydrogen_constraints.extend(constraints)

    def _parse_hydrogen_constraint(self, constraint):
        symbol = constraint.pop("symbol", None)

        if symbol is None:
            raise RuntimeError("Invalid hydrogen constraint: " +
                               str(constraint) +
                               " - No symbol provided")
        return {symbol: constraint}

    def to_dict(self):
        return {"BeamRadius": self._beam_radius,
                "SampleHeight": self._sample_height,
                "SampleWidth": self._sample_width,
                "SampleDepth": self._sample_depth,
                "SampleDensity": self._sample_density,
                "Seed": self._seed,
                "NumScatters": self._num_scatters,
                "NumRuns": self._num_runs,
                "NumEvents": self._num_events,
                "SmoothNeighbours": self._smooth_neighbours,
                "HydrogenConstraints": self._hydrogen_constraints.to_dict()}


# -----------------------------------------------------------------------------------------

class VesuvioTOFFitHelper(object):
    """
    A helper class for executing a VesuvioTOFFit.

    Attributes:
        _background              The background to use in the fit
        _intensity_constraints   The intensity constraints to use in the fit
        _ties                    The ties to use in the fit
        _max_iterations          The max number of iterations in which to attempt to find
                                 a peak using the minimizer
        _minimizer               The minimizer to use in the fit
    """
    def __init__(self, background, intensity_constraints, ties, max_iterations, minimizer):
        self._background = background
        self._intensity_constraints = intensity_constraints
        self._ties = ties
        self._max_iterations = max_iterations
        self._minimizer = minimizer

    def __call__(self, **fit_args):
        return VesuvioTOFFit(Background=self._background,
                             IntensityConstraints=self._intensity_constraints,
                             Ties=self._ties,
                             MaxIterations=self._max_iterations,
                             Minimizer=self._minimizer,
                             **fit_args)


# -----------------------------------------------------------------------------------------

class VesuvioCorrectionsHelper(object):
    """
    A helper class for executing VesuvioCorrections.

    Attributes:
        _gamma_correct          Boolean specifying whether to calculate gamma
                                background corrections
        _multiple_scattering    Boolean specifying whether to calculate multiple
                                scattering corrections
        _gamma_background_scale The value to scale the gamma corrections by
        _container_scale        The value to scale the container corrections by
        _intensity_constraints  The intensity constraints to use in the corrections
                                calculation
    """
    def __init__(self, gamma_correct, multiple_scattering, gamma_background_scale,
                 container_scale, intensity_constraints):
        self._gamma_correct = gamma_correct
        self._multiple_scattering = multiple_scattering
        self._gamma_background_scale = gamma_background_scale
        self._container_scale = container_scale
        self._intensity_constraints = intensity_constraints

    def __call__(self, **correction_args):
        return VesuvioCorrections(GammaBackground=self._gamma_correct,
                                  IntensityConstraints=self._intensity_constraints,
                                  MultipleScattering=self._multiple_scattering,
                                  GammaBackgroundScale=self._gamma_background_scale,
                                  ContainerScale=self._container_scale,
                                  **correction_args)


# -----------------------------------------------------------------------------------------

class VesuvioConstraints(object):
    """
    A class for parsing and storing a set of constraints for the Vesuvio Fit Routine

    Attributes:
        _name           The name of this set of constraints
        _parser         The parser to be used
        _constraints    The set of parsed constraints
    """
    def __init__(self, name, parser):
        self._name = name
        self._parser = parser
        self._constraints = dict()

    def __iter__(self):
        return self._constraints.__iter__()

    def __contains__(self, key):
        return self._constraints.__contains__(key)

    def add(self, key, constraint):
        self._constraints[key] = constraint

    def extend(self, constraints):
        try:
            if not isinstance(constraints, dict):
                for constraint in constraints:
                    self._constraints.update(self._parser(constraint))
            else:
                self._constraints.update(self._parser(constraints))
        except AttributeError:
            raise RuntimeError(self._name + ": Constraints are incorrectly formatted.")

    def to_dict(self):
        return dict(self._constraints)


# -----------------------------------------------------------------------------------------

class MassProfileCollection(object):
    """
    A class for storing a collection of mass values and profiles.

    Attributes:
        _material_builder       A material builder for retrieving mass
                                values from symbols
        _profiles               A list containing the mass profiles for
                                each spectrum index
        _index_to_symbol_map    A map from the index of a mass in the
                                mass values list to its chemical symbol,
                                where specified
        _mass_values            A list of the mass values
        _profile_functions      The function strings corresponding to the
                                mass  profiles
    """
    def __init__(self, profiles):
        self._material_builder = MaterialBuilder()
        self._profiles = [profiles]
        self._index_to_symbol_map = {}
        self._set_profile_functions(self._profiles)
        self._mass_values = []
        self._extract_mass_values(profiles, True)

    def profile_functions(self, index=0, symbol_filter=None):
        profile_functions = self._profile_functions[min(index, len(self._profile_functions)-1)]

        if symbol_filter is None:
            return ';'.join(profile_functions)
        else:
            return ';'.join(self._apply_symbol_filter(profile_functions, symbol_filter))

    def masses(self, symbol_filter=None):
        if symbol_filter is None:
            return self._mass_values
        else:
            return self._apply_symbol_filter(self._mass_values, symbol_filter)

    def index_to_symbol_map(self):
        return self._index_to_symbol_map

    def update_profiles_from_workspace(self, params_ws, symbol_filter=None):
        function_regex = re.compile("f([0-9]+).([A-z0-9_]+)")
        param_labels = params_ws.getAxis(1).extractValues()
        num_masses = len(self._mass_values)
        self._resize_profiles(params_ws.blocksize())

        if symbol_filter is None:
            mass_indices = range(num_masses)
        else:
            mass_indices = self._apply_symbol_filter(self._mass_values, symbol_filter,
                                                     create_index_list=True)

        for idx in range(params_ws.blocksize()):

            for param_idx, param in enumerate(param_labels):
                if param != "Cost function value":
                    param_re = function_regex.match(param)
                    mass_idx = mass_indices[int(param_re.group(1))]

                    if mass_idx >= num_masses:
                        continue

                    param_name = param_re.group(2).lower()
                    self._update_profile_parameter(param_name, params_ws.dataY(param_idx)[idx],
                                                   idx, mass_idx)
        self._set_profile_functions(self._profiles)

    def _apply_symbol_filter(self, collection, symbol_filter, create_index_list=False):
        filtered = []
        for idx, item in enumerate(collection):
            if idx not in self._index_to_symbol_map or \
                    symbol_filter(self._index_to_symbol_map[idx]):
                filtered.append(idx if create_index_list else item)
        return filtered

    def _set_profile_functions(self, profiles_list):
        self._profile_functions = [[self._extract_function(profile) for profile in profiles]
                                   for profiles in profiles_list]

    def _extract_function(self, profile_properties):
        function_props = []

        try:
            function_type = profile_properties['function']
            function_props.append(self._create_function_parameter('function', function_type))
        except KeyError:
            raise RuntimeError("Function not specified in profile: "
                               + str(profile_properties))

        for key, value in profile_properties.items():
            if key != 'function':
                function_props.append(self._create_function_parameter(key, value))
        return ','.join(function_props)

    def _create_function_parameter(self, param_key, param_value):
        if param_key != 'symbol':
            return "{0}={1}".format(param_key, param_value)
        else:
            return "{0}={1}".format('value', self._mass_from_chemical_symbol(param_value))

    def _extract_mass_values(self, profiles, create_index_to_symbol_map=False):
        for profile in profiles:
            mass_value = self._extract_mass_value(profile)
            self._mass_values.append(mass_value)
            if create_index_to_symbol_map:
                self._index_to_symbol_map[len(self._mass_values)] = mass_value

    def _extract_mass_value(self, profile_properties):
        mass_value = profile_properties.get('value', None)
        if mass_value is None:
            symbol = profile_properties.get('symbol', None)

            if symbol is None:
                raise RuntimeError('Invalid mass specified - ' + str(profile_properties)
                                   + " - either 'value' or 'symbol' must be given.")

            try:
                mass_value = self._mass_from_chemical_symbol(symbol)
            except BaseException as exc:
                raise RuntimeError('Error when parsing mass - ' + str(profile_properties) + ": "
                                   + "\n" + str(exc))
        return mass_value

    def _resize_profiles(self, size):
        while len(self._profiles) <= size:
            profiles = []
            for idx in range(len(self._mass_values)):
                profiles.append({'function' : self._profiles[0][idx]['function']})
            self._profiles.append(profiles)

    def _update_profile_parameter(self, param_name, param_value, index, mass_index):
        if param_name == 'width' and \
                isinstance(self._profiles[index][mass_index].get(param_name, None), list):
            self._profiles[index][mass_index][param_name][1] = param_value
        else:
            self._profiles[index][mass_index][param_name] = param_value

    def _mass_from_chemical_symbol(self, chemical_symbol):
        return self._material_builder.setFormula(chemical_symbol).build().relativeMolecularMass()


# -----------------------------------------------------------------------------------------

class VesuvioTOFFitOutput(object):

    def __init__(self, vesuvio_input, verbose_corrections):
        self._verbose_corrections = verbose_corrections
        self._input = vesuvio_input
        self._group_name = vesuvio_input.sample_runs + "_result"
        self._num_spectra = vesuvio_input.sample_data.getNumberHistograms()
        self._result_ws = None
        self._pre_correct_pars_workspace = None
        self._pars_workspace = None
        self._fit_workspace = None
        self._iteration = 1
        self._result_workspaces = []
        self._data_workspaces = []
        self._chi2_values = []
        self._previous_chi2_values = []
        self._fit_ws_name = ""

    def contains_fit(self):
        return bool(self._chi2_values)

    def chi2_values(self):
        return self._chi2_values

    def prefit_parameters_workspace(self):
        return self._pre_correct_pars_workspace

    def parameters_workspace(self):
        return self._pars_workspace

    def result_workspace(self):
        if self._result_ws is None:
            self._result_ws = self._create_result_workspace()
        return self._result_ws

    def change_in_cost_function(self):
        last_chi2 = np.array(self._previous_chi2_values)
        chi2 = np.array(self._chi2_values)
        chi2_delta = last_chi2 - chi2
        return np.abs(np.max(chi2_delta))

    def set_fit_workspace_name(self, fit_ws_name):
        self._fit_ws_name = fit_ws_name

    def add_chi2_value(self, chi2_value):
        self._chi2_values.append(chi2_value)

    def update_prefit_parameters_workspace(self, table_workspace, workspace_index):
        self._pre_correct_pars_workspace = self._update_workspace(self._pre_correct_pars_workspace,
                                                                  table_workspace, workspace_index)

    def update_parameters_workspace(self, table_workspace, workspace_index):
        self._pars_workspace = self._update_workspace(self._pars_workspace, table_workspace, workspace_index)

    def update_fit_workspace(self, table_workspace, workspace_index):
        self._fit_workspace = self._update_workspace(self._fit_workspace, table_workspace, workspace_index)

    def increment_iteration(self):
        self._result_ws = None
        self._pre_correct_pars_workspace = None
        self._pars_workspace = None
        self._fit_workspace = None
        self._result_workspaces = []
        self._data_workspaces = []
        self._previous_chi2_values = self._chi2_values
        self._chi2_values = []
        self._iteration += 1

    def add_output_to_ads(self, workspace_index, correction_workspaces, corrected_workspaces):

        if self._verbose_corrections:
            self._create_verbose_output(workspace_index, correction_workspaces, corrected_workspaces)

        sample_runs = self._input.sample_runs
        params_pre_corr = sample_runs + "_params_pre_correction_iteration_" + str(self._iteration)
        params_name = sample_runs + "_params_iteration_" + str(self._iteration)
        fit_name = sample_runs + "_correction_fit_scale_iteration_" + str(self._iteration)
        mtd.addOrReplace(params_pre_corr, self._pre_correct_pars_workspace)
        mtd.addOrReplace(params_name, self._pars_workspace)
        mtd.addOrReplace(fit_name, self._fit_workspace)

    def _update_workspace(self, workspace, table_workspace, workspace_index):
        if workspace is None:
            workspace = self._create_param_workspace(self._num_spectra, table_workspace)

        sample_data = self._input.sample_data
        spectrum = sample_data.getSpectrum(workspace_index).getSpectrumNo()
        current_spectrum = 'spectrum_' + str(spectrum)
        self._update_fit_params(workspace, workspace_index, table_workspace, current_spectrum)
        return workspace

    def _create_verbose_output(self, workspace_index, correction_workspaces, corrected_workspaces):
        output_workspaces = correction_workspaces.getNames()
        output_workspaces += corrected_workspaces.getNames()
        UnGroupWorkspace(correction_workspaces)
        UnGroupWorkspace(corrected_workspaces)

        for workspace in output_workspaces:
            self._group_name = self._input.sample_runs + '_iteration' + str(self._iteration)
            name = self._group_name + '_' + workspace.split('_')[1] + '_' + workspace.split('_')[-1]
            self._result_workspaces.append(name)

            if workspace_index == 0:
                RenameWorkspace(InputWorkspace=workspace, OutputWorkspace=name)
            else:
                ConjoinWorkspaces(InputWorkspace1=name, InputWorkspace2=workspace)

    def _create_result_workspace(self):
        output_groups = []

        if self._result_workspaces:
            output_groups.append(GroupWorkspaces(InputWorkspaces=self._result_workspaces,
                                                 OutputWorkspace=self._group_name))
        if self._data_workspaces:
            output_groups.append(GroupWorkspaces(InputWorkspaces=self._data_workspaces,
                                                 OutputWorkspace=self._group_name + '_data'))
        else:
            output_groups.append(self._fit_ws_name)

        if len(output_groups) > 1:
            return output_groups
        else:
            return output_groups[0]

    def _update_fit_params(self, params_ws, spec_idx, params_table, name):
        params_ws.getAxis(0).setLabel(spec_idx, name)
        for idx in range(params_table.rowCount()):
            params_ws.dataX(idx)[spec_idx] = spec_idx
            params_ws.dataY(idx)[spec_idx] = params_table.column('Value')[idx]
            params_ws.dataE(idx)[spec_idx] = params_table.column('Error')[idx]

    def _create_param_workspace(self, num_spec, param_table):
        num_params = param_table.rowCount()
        param_workspace = WorkspaceFactory.Instance().create("Workspace2D",
                                                             num_params, num_spec,
                                                             num_spec)
        x_axis = TextAxis.create(num_spec)
        param_workspace.replaceAxis(0, x_axis)

        vert_axis = TextAxis.create(num_params)
        for idx, param_name in enumerate(param_table.column('Name')):
            vert_axis.setLabel(idx, param_name)
        param_workspace.replaceAxis(1, vert_axis)

        return param_workspace


# --------------------------------------------------------------------------------
# Private Functions
# --------------------------------------------------------------------------------

def _create_background_str(background_flags):
    """
    Create a string suitable for the algorithms out of the background flags
    :param background_flags: A dict for the background (can be None)
    :return: A string to pass to the algorithm
    """
    if background_flags:
        background_props = ["function={0}".format(background_flags["function"])]
        del background_flags["function"]
        for key, value in iteritems(background_flags):
            background_props.append("{0}={1}".format(key, value))
        background_str = ",".join(background_props)
    else:
        background_str = ""

    return background_str


def _create_intensity_constraint_str(intensity_constraints):
    """
    Create a string suitable for the algorithms out of the intensity constraint flags
    :param intensity_constraints: A list of lists for the constraints (can be None)
    :return: A string to pass to the algorithm
    """
    if intensity_constraints:
        if not isinstance(intensity_constraints[0], list):
            intensity_constraints = [intensity_constraints]
        # Make each element a string and then join them together
        intensity_constraints = [str(c) for c in intensity_constraints]
        intensity_constraints_str = ";".join(intensity_constraints)
    else:
        intensity_constraints_str = ""

    return intensity_constraints_str


def _create_user_defined_ties_str(masses):
    """
    Creates the internal ties for each mass profile as defined by the user to be used when fitting the data
    @param masses   :: The mass profiles for the data which contain the the ties
    @return         :: A string to be passed as the Ties input to fitting
    """
    user_defined_ties = []
    for index, mass in enumerate(masses):
        if 'ties' in mass:
            ties = mass['ties'].split(',')
            function_identifier = 'f' + str(index) + '.'
            for t in ties:
                tie_str = function_identifier + t
                equal_pos = tie_str.index('=') + 1
                tie_str = tie_str[:equal_pos] + function_identifier + tie_str[equal_pos:]
                user_defined_ties.append(tie_str)
    user_defined_ties = ','.join(user_defined_ties)
    return user_defined_ties


def ignore_hydrogen_filter(symbol):
    """
    A filter which returns False if the specified chemical symbol is
    for hydrogen and True otherwise.

    :param symbol:  The chemical symbol.
    :return:        True if the specified chemical symbol is for
                    hydrogen, False otherwise.
    """
    return symbol != 'H'
