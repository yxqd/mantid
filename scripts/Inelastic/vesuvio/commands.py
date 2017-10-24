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
from mantid.api import (AnalysisDataService, MatrixWorkspace, WorkspaceFactory, TextAxis)
from mantid.kernel import MaterialBuilder
from vesuvio.instrument import VESUVIO

import mantid.simpleapi as ms


# --------------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------------

def fit_tof(runs, flags, iterations=1, convergence_threshold=None):
    vesuvio_loader = VesuvioLoadHelper(flags['diff_mode'], flags['fit_mode'],
                                       flags['ip_file'], flags['bin_parameters'])
    vesuvio_input = VesuvioTOFFitInput(flags['runs'], flags['container_runs'],
                                       flags['spectra'], vesuvio_loader)
    ms_helper = None
    if flags.get('ms_enabled', False):
        ms_helper = VesuvioMSHelper(**flags['ms_flags'])

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
    return vesuvio_fit_routine(vesuvio_input, iterations, convergence_threshold,
                               flags.get('output_verbose_corrections', False))

# -----------------------------------------------------------------------------------------

class VesuvioTOFFitInput(object):

    def __init__(self, sample_runs, container_runs, spectra, loader):
        self.sample_runs = sample_runs
        self.container_runs = container_runs
        self.sample_data = self._load_data(sample_runs, spectra, loader)
        self.container_data = self._load_data(container_runs, spectra, loader)
        self.spectra = spectra
        self._back_scattering = self._is_back_scattering_spectra(spectra)
        mtd.addOrReplace(self._tof_workspace_suffix(sample_runs, spectra),
                         self.sample_data)
        mtd.addOrReplace(self._tof_workspace_suffix(container_runs, spectra),
                         self.container_data)

    def using_back_scattering_spectra(self):
        return self._back_scattering

    def _is_back_scattering_spectra(self, spectra):
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

    def _load_data(self, runs, spectra, loader):
        if isinstance(runs, MatrixWorkspace):
            return runs
        else:
            return loader(runs, spectra)

    def _tof_workspace_suffix(self, runs, spectra):
        return runs + "_" + spectra + "_tof"

# -----------------------------------------------------------------------------------------

class VesuvioTOFFitRoutine(object):

    def __init__(self, ms_helper, fit_helper, corrections_helper, mass_profile_collection, fit_mode):
        self._ms_helper = ms_helper
        self._fit_helper = fit_helper
        self._corrections_helper = corrections_helper
        self._mass_profile_collection = mass_profile_collection
        self._fit_mode = fit_mode

    def __call__(self, vesuvio_input, iterations, convergence_threshold, verbose_output=False):
        if iterations < 1:
            return ValueError('Must perform at least one iteration')

        tof_iteration = VesuvioTOFFitRoutineIteration(self._ms_helper, self._fit_helper,
                                                      self._mass_profile_collection)

        back_scattering = vesuvio_input.using_back_scattering()
        ignore_hydrogen_filter = lambda x : x is not 'H'
        update_filter = ignore_hydrogen_filter if back_scattering else None
        exit_iteration = 0
        previous_results = None

        for iteration in range(1, iterations+1):
            if previous_results is not None:
                self._mass_profile_collection.update_profiles_from_workspace(previous_results[2],
                                                                             update_filter)
            print("=== Iteration {0} out of a possible {1}".format(iteration, iterations))
            results = tof_iteration(vesuvio_input, iteration, verbose_output)
            exit_iteration += 1

            if previous_results is not None and convergence_threshold is not None:
                cost_function_change = self._calculate_change_in_cost_function(results,
                                                                               previous_results)
                print("Cost function change: {0}".format(cost_function_change))

                if cost_function_change <= convergence_threshold:
                    print("Stopped at iteration {0} due to minimal change in cost function"
                          .format(exit_iteration))
                    return results[0], results[2], results[3], exit_iteration
            previous_results = results
        return previous_results[0], previous_results[2], previous_results[3], exit_iteration

    def _calculate_change_in_cost_function(self, results, previous_results):
        last_chi2 = np.array(previous_results[3])
        chi2 = np.array(results[3])
        chi2_delta = last_chi2 - chi2
        return np.abs(np.max(chi2_delta))

# ------------------------------------------------------------------------------------------------------

class VesuvioTOFFitRoutineIteration(object):

    def __init__(self, ms_helper, fit_helper, corrections_helper, mass_profile_collection, fit_mode):
        self._fit_mode = fit_mode
        self._ms_helper = ms_helper
        self._fit_helper = fit_helper
        self._corrections_helper = corrections_helper
        self._mass_profile_collection = mass_profile_collection

    def __call__(self, vesuvio_input, iteration, verbose_output=False):
        sample_runs = str(vesuvio_input.sample_runs)
        sample_data = vesuvio_input.sample_data
        container_data = vesuvio_input.container_data
        num_spectra = sample_data.getNumberHistograms()

        ignore_hydrogen_filter = lambda x : x is not 'H'
        fit_filter = ignore_hydrogen_filter if vesuvio_input.using_back_scattering_spectra() else None
        all_mass_values = self._mass_profile_collection.masses()
        fit_mass_values = self._mass_profile_collection.masses(fit_filter)

        # Outputs
        group_name = sample_runs + '_result'
        fit_ws_name = ""
        pre_correct_pars_workspace = None
        pars_workspace = None
        fit_workspace = None
        output_groups = []
        chi2_values = []
        data_workspaces = []
        result_workspaces = []

        for index in range(num_spectra):
            all_profiles = self._mass_profile_collection.profiles(index)
            fit_profiles = self._mass_profile_collection.profiles(index, fit_filter)

            # Calculate pre-fit to retrieve parameter approximations for corrections
            suffix = self._create_fit_workspace_suffix(index, sample_data, self._fit_mode, iteration)
            fit_ws_name = "__vesuvio_corrections_fit"
            pre_params_ws_name = sample_runs + "_params_pre_correction" + suffix
            self._fit_helper(sample_data, index, fit_mass_values, fit_profiles,
                             fit_ws_name, pre_params_ws_name)
            ms.DeleteWorkspace(fit_ws_name)

            # Calculate and apply vesuvio corrections
            linear_corrections_fit_params_name = sample_runs + "_correction_fit_scale" + suffix
            corrected_data_name = sample_runs + "_tof_corrected" + suffix
            corrections_args = {'FitParameters' : pre_params_ws_name}
            corrections_args.update(self._ms_helper.to_dict())
            if verbose_output:
                corrections_args['CorrectionsWorkspaces'] = sample_runs + "_correction" + suffix
                corrections_args['CorrectedWorkspaces'] = sample_runs + "_corrected" + suffix
            if container_data is not None:
                corrections_args['ContainerWorkspace'] = container_data
            self._corrections_helper(sample_data, index, all_mass_values, all_profiles,
                                     self._mass_profile_collection.index_to_symbol_map(),
                                     corrected_data_name, linear_corrections_fit_params_name,
                                     **corrections_args)

            # Calculate final fit
            fit_ws_name = sample_runs + "_data" + suffix
            params_ws_name = sample_runs + "_params" + suffix
            fit_result = self._fit_helper(corrected_data_name, 0, fit_mass_values,
                                          fit_profiles, fit_ws_name, params_ws_name)
            chi2_values.append(fit_result[-1])
            ms.DeleteWorkspace(corrected_data_name)

            # Post-processing of parameter tables - update workspaces using param tables
            if pre_correct_pars_workspace is None:
                pre_correct_pars_workspace = \
                    self._create_param_workspace(num_spectra, mtd[pre_params_ws_name])
            if pars_workspace is None:
                pars_workspace = \
                    self._create_param_workspace(num_spectra, mtd[params_ws_name])
            if fit_workspace is None:
                fit_workspace = \
                    self._create_param_workspace(num_spectra, mtd[linear_corrections_fit_params_name])

            spectrum = sample_data.getSpectrum(index).getSpectrumNo()
            current_spectrum = 'spectrum_' + str(spectrum)
            self._update_fit_params(pre_correct_pars_workspace, index, mtd[pre_params_ws_name],
                                    current_spectrum)
            self._update_fit_params(pars_workspace, index, mtd[params_ws_name], current_spectrum)
            self._update_fit_params(fit_workspace, index, mtd[linear_corrections_fit_params_name],
                                    current_spectrum)

            # Delete parameter tables after use
            ms.DeleteWorkspace(pre_params_ws_name)
            ms.DeleteWorkspace(params_ws_name)
            ms.DeleteWorkspace(linear_corrections_fit_params_name)

            # Process spectrum group
            # Note the ordering of operations here gives the order in the WorkspaceGroup
            output_workspaces = []
            data_workspaces.append(fit_ws_name)
            if verbose_output:
                output_workspaces += mtd[corrections_args["CorrectionWorkspaces"]].getNames()
                output_workspaces += mtd[corrections_args["CorrectedWorkspaces"]].getNames()
                ms.UnGroupWorkspace(corrections_args["CorrectionWorkspaces"])
                ms.UnGroupWorkspace(corrections_args["CorrectedWorkspaces"])

                for workspace in output_workspaces:
                    group_name = sample_runs + '_iteration_' + str(iteration)
                    name = group_name + '_' + workspace.split('_')[1] + '_' + workspace.split('_')[-1]
                    result_workspaces.append(name)
                    if index == 0:
                        ms.RenameWorkspace(InputWorkspace=workspace, OutputWorkspace=name)
                    else:
                        ms.ConjoinWorkspaces(InputWorkspace1=name, InputWorkspace2=workspace)

            # Output the parameter workspaces
            params_pre_corr = sample_runs + "_params_pre_correction_iteration_" + str(iteration)
            params_name = sample_runs + "_params_iteration_" + str(iteration)
            fit_name = sample_runs + "_correction_fit_scale_iteration_" + str(iteration)
            mtd.addOrReplace(params_pre_corr, pre_correct_pars_workspace)
            mtd.addOrReplace(params_name, pars_workspace)
            mtd.addOrReplace(fit_name, fit_workspace)

        if result_workspaces:
            output_groups.append(ms.GroupWorkspaces(InputWorkspaces=result_workspaces,
                                                    OutputWorkspace=group_name))

        if data_workspaces:
            output_groups.append(ms.GroupWorkspaces(InputWorkspaces=data_workspaces,
                                                    OutputWorkspace=group_name + '_data'))
        else:
            output_groups.append(fit_ws_name)

        if len(output_groups) > 1:
            result_ws = output_groups
        else:
            result_ws = output_groups[0]
        return result_ws, pre_correct_pars_workspace, pars_workspace, chi2_values

    def _create_fit_workspace_suffix(self, index, tof_data, fit_mode, spectra, iteration=None):
        if fit_mode == "bank":
            suffix = "_" + spectra + "_bank_" + str(index + 1)
        else:
            spectrum = tof_data.getSpectrum(index)
            suffix = "_spectrum_" + str(spectrum.getSpectrumNo())

        if iteration is not None:
            suffix += "_iteration_" + str(iteration)

        return suffix

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

    def _update_fit_params(self, params_ws, spec_idx, params_table, name):
        params_ws.getAxis(0).setLabel(spec_idx, name)
        for idx in range(params_table.rowCount()):
            params_ws.dataX(idx)[spec_idx] = spec_idx
            params_ws.dataY(idx)[spec_idx] = params_table.column('Value')[idx]
            params_ws.dataE(idx)[spec_idx] = params_table.column('Error')[idx]

# ------------------------------------------------------------------------------------------------------

class VesuvioLoadHelper(object):
    """
    A helper class for loading Vesuvio data from the input of a user script.
    """

    def __init__(self, diff_mode, fit_mode, param_file, rebin_params=None):
        self._diff_mode = diff_mode
        self._fit_mode = fit_mode
        self._param_file = param_file
        self._rebin_params = rebin_params
        self._instrument = VESUVIO()

    def __call__(self, runs, spectra):
        return self.load_and_crop_runs(runs, spectra)

    def load_and_crop_runs(self, runs, spectra):
        sum_spectra = (self._fit_mode == 'bank')
        loaded = self.load_runs(runs, spectra, sum_spectra)
        cropped = self._crop_workspace(loaded)

        if self._rebin_params is not None:
            return self._rebin_workspace(cropped)
        else:
            return cropped

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
        load_alg = AlgorithmManager.createUnmanaged("LoadVesuvio")
        load_alg.setProperty("Filename", runs)
        load_alg.setProperty("Mode", self._diff_mode)
        load_alg.setProperty("InstrumentParFile", self._param_file)
        load_alg.setProperty("SpectrumList", spectra)
        load_alg.setProperty("SumSpectra", sum_spectra)
        load_alg.setProperty("OutputWorkspace", "loaded")
        load_alg.execute()
        return load_alg.getProperty("OutputWorkspace").value

    def _crop_workspace(self, workspace):
        crop_alg = AlgorithmManager.createUnmanaged("CropWorkspace")
        crop_alg.setProperty("InputWorkspace", workspace)
        crop_alg.setProperty("XMin", self._instrument.tof_range[0])
        crop_alg.setProperty("XMax", self._instrument.tof_range[1])
        crop_alg.setProperty("OutputWorkspace", "cropped")
        crop_alg.execute()
        return crop_alg.getProperty("OutputWorkspace").value

    def _rebin_workspace(self, workspace):
        rebin_alg = AlgorithmManager.createUnmanaged("Rebin")
        rebin_alg.setProperty("InputWorkspace", workspace)
        rebin_alg.setProperty("Params", self._rebin_params)
        rebin_alg.setProperty("OutputWorkspace", "rebinned")
        rebin_alg.execute()
        return rebin_alg.getProperty("OutputWorkspace")

# ------------------------------------------------------------------------------------------------------

class VesuvioMSHelper(object):

    def __init__(self, BeamRadius=2.5, SampleHeight=5.0, SampleWidth=5.0, SampleDepth=5.0,
                 SampleDensity=1.0, Seed=123456789, NumScatters=3, NumRuns=10, NumEvents=50000,
                 SmoothNeighbours=3):
        self.beam_radius = BeamRadius
        self.sample_height = SampleHeight
        self.sample_width = SampleWidth
        self.sample_depth = SampleDepth
        self.sample_density = SampleDensity
        self.seed = Seed
        self.num_scatters = NumScatters
        self.num_runs = NumRuns
        self.num_events = NumEvents
        self.smooth_neighbours = SmoothNeighbours
        parser = lambda x : self._parse_hydrogen_constraint(x)
        self._hydrogen_constraints = VesuvioConstraints("HydrogenConstraints", parser)

    def add_hydrogen_constraints(self, constraints):
        self._hydrogen_constraints.extend(constraints)

    def _parse_hydrogen_constraint(self, constraint):
        symbol = constraint.pop("symbol", None)

        if symbol is None:
            raise RuntimeError("Invalid hydrogen constraint: " +
                               str(constraint) +
                               " - No symbol provided")
        return {symbol : constraint}

    def to_dict(self):
        return {"BeamRadius" : self.beam_radius,
                "SampleHeight" : self.sample_height,
                "SampleWidth" : self.sample_width,
                "SampleDepth" : self.sample_depth,
                "SampleDensity" : self.sample_density,
                "Seed" : self.seed,
                "NumScatters" : self.num_scatters,
                "NumRuns" : self.num_runs,
                "NumEvents" : self.num_events,
                "SmoothNeighbours" : self.smooth_neighbours,
                "HydrogenConstraints" : self._hydrogen_constraints.to_dict()}

# -----------------------------------------------------------------------------------------

class VesuvioTOFFitHelper(object):

    def __init__(self, background, intensity_constraints, ties, max_iterations, minimizer):
        self._background = background
        self._intensity_constraints = intensity_constraints
        self._ties = ties
        self._max_iterations = max_iterations
        self._minimizer = minimizer

    def __call__(self, input_workspace, workspace_index, masses, mass_profiles, output_ws_name, params_ws_name):
        return ms.VesuvioTOFFit(InputWorkspace=input_workspace,
                                WorkspaceIndex=workspace_index,
                                Masses=masses,
                                MassProfiles=mass_profiles,
                                Background=self._background,
                                IntensityConstraints=self._intensity_constraints,
                                Ties=self._ties,
                                OutputWorkspace=output_ws_name,
                                FitParameters=params_ws_name,
                                MaxIterations=self._max_iterations,
                                Minimizer=self._minimizer)

# -----------------------------------------------------------------------------------------

class VesuvioCorrectionsHelper(object):

    def __init__(self, gamma_correct, multiple_scattering, gamma_background_scale,
                 container_scale, intensity_constraints):
        self._gamma_correct = gamma_correct
        self._multiple_scattering = multiple_scattering
        self._gamma_background_scale = gamma_background_scale
        self._container_scale = container_scale
        self._intensity_constraints = intensity_constraints

    def __call__(self, input_workspace, workspace_index, masses, mass_profiles,
                 index_to_symbol_map, output_workspace, linear_fit_results,
                 **corrections_args):
        return ms.VesuvioCorrections(InputWorkspace=input_workspace,
                                     OutputWorkspace=output_workspace,
                                     LinearFitResult=linear_fit_results,
                                     WorkspaceIndex=workspace_index,
                                     GammaBackground=self._gamma_correct,
                                     Masses=masses,
                                     MassProfiles=mass_profiles,
                                     MassIndexToSymbolMap=index_to_symbol_map,
                                     IntensityConstraints=self._intensity_constraints,
                                     MultipleScattering=self._multiple_scattering,
                                     GammaBackgroundScale=self._gamma_background_scale,
                                     ContainerScale=self._container_scale,
                                     **corrections_args)

# -----------------------------------------------------------------------------------------

class VesuvioConstraints(object):

    def __init__(self, name, parser):
        self._name = name
        self._parser = parser
        self._constraints = dict()

    def __iter__(self):
        return self._constraints.__iter__()

    def __contains__(self, key):
        return self._constraints.__contains__(key)

    def key_complement(self, keys):
        return [key for key in self._constraints if key not in keys]

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

    def __init__(self, profiles):
        self._profiles = [profiles]
        self._index_to_symbol_map = {}
        self._set_profile_functions(self._profiles)
        self._mass_values = [self._extract_mass_value(profile, True) for profile in profiles]
        self._material_builder = MaterialBuilder()

    def profile_functions(self, index=0, symbol_filter=None):
        if symbol_filter is None:
            return ';'.join(self._profile_functions[index])
        else:
            return ';'.join(self._apply_symbol_filter(self._profile_functions[index], symbol_filter))

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
        mass_indices = self._apply_symbol_filter(self._mass_values[0], symbol_filter,
                                                 create_index_list=True)
        num_masses = len(self._mass_values[0])
        self._resize_profiles(params_ws.blocksize(), num_masses)

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

        self._profile_functions = [[self._extract_function(profile) for profile in profiles]
                                   for profiles in self._profiles]

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
        function_props = [self._extract_function_parameter(key, value) for key, value in
                          profile_properties.items()]
        return ','.join(function_props)

    def _extract_function_parameter(self, param_key, param_value):
        if param_key != 'symbol':
            return "{0}={1}".format(param_key, param_value)
        else:
            return "{0}={1}".format('value', self._mass_from_chemical_symbol(param_key))

    def _extract_mass_value(self, profile_properties, create_index_to_symbol_map=False):
        mass_value = profile_properties.get('value', None)
        if mass_value is None:
            symbol = profile_properties.get('symbol', None)

            if symbol is None:
                raise RuntimeError('Invalid mass specified - ' + str(profile_properties)
                                   + " - either 'value' or 'symbol' must be given.")

            try:
                mass_value = self._mass_from_chemical_symbol(symbol)
                if create_index_to_symbol_map:
                    self._index_to_symbol_map[len(self._mass_values)] = mass_value
            except BaseException as exc:
                raise RuntimeError('Error when parsing mass - ' + str(profile_properties) + ": "
                                   + "\n" + str(exc))
        return mass_value

    def _resize_profiles(self, dim1, dim2):
        while len(self._profiles) <= dim1.blocksize():
            profiles = []
            while len(profiles) <= dim2:
                profiles.append({})
            self._profiles.append(profiles)

    def _update_profile_parameter(self, param_name, param_value, index, mass_index):
        if param_name == 'width' and \
                isinstance(self._profiles[index][mass_index].get(param_name, None), list):
            self._profiles[index][mass_index][param_name][1] = param_value
        else:
            self._profiles[index][mass_index][param_name] = param_value

    def _mass_from_chemical_symbol(self, chemical_symbol):
        return self._material_builder.setFormula(chemical_symbol).build().relativeMolecularMass()

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
