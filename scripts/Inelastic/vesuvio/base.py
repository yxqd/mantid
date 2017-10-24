# pylint: disable=too-few-public-methods,redefined-builtin
from __future__ import (absolute_import, division, print_function)
from six import iteritems

from mantid.api import Algorithm, AlgorithmManager
from vesuvio.instrument import VESUVIO


class VesuvioBase(Algorithm):

    # There seems to be a problem with Python algorithms
    # defining a __init__ method
    _INST = None

    def _execute_child_alg(self, name, **kwargs):
        """Execute an algorithms as a child.
        By default all outputs are returned but this can be limited
        by providing the return_values=[] keyword
        """
        alg = self.createChildAlgorithm(name)
        # For Fit algorithm, Function & InputWorkspace have to
        # be set first and in that order.
        if name == 'Fit':
            for key in ('Function', 'InputWorkspace'):
                alg.setProperty(key, kwargs[key])
                del kwargs[key]

        ret_props = None
        if 'return_values' in kwargs:
            ret_props = kwargs['return_values']
            if type(ret_props) is str:
                ret_props = [ret_props]
            del kwargs['return_values']

        for name, value in iteritems(kwargs):
            alg.setProperty(name, value)
        alg.execute()

        # Assemble return values
        if ret_props is None:
            # This must be AFTER execute just in case that attached more
            # output properties
            ret_props = alg.outputProperties()
        outputs = []
        for name in ret_props:
            outputs.append(alg.getProperty(name).value)
        if len(outputs) == 1:
            return outputs[0]
        else:
            return tuple(outputs)

# -----------------------------------------------------------------------------------------
# Helper to translate from an table workspace to a dictionary. Should be on the workspace
# really ...
# -----------------------------------------------------------------------------------------


class TableWorkspaceDictionaryFacade(object):
    """
    Allows an underlying table workspace to be treated like a read-only dictionary
    """

    def __init__(self, held_object):
        self._table_ws = held_object

    def __getitem__(self, item):
        for row in self._table_ws:
            if row['Name'] == item:
                return row['Value']
        #endfor
        raise KeyError(str(item))

    def __contains__(self, item):
        for row in self._table_ws:
            if row['Name'] == item:
                return True

        return False

# -----------------------------------------------------------------------------------------

class MassProfileCollection(object):

    def __init__(self):
        self._profiles = []
        self._mass_values = []
        self._index_to_symbol_map = {}

    def profiles(self, index=0, symbol_filter=None):
        if symbol_filter is None:
            return ';'.join(self._profiles[index])
        else:
            return ';'.join(self._apply_symbol_filter(self._profiles[index], symbol_filter))

    def masses(self, index=0, symbol_filter=None):
        if symbol_filter is None:
            return self._mass_values[index]
        else:
            return self._apply_symbol_filter(self._mass_values[index], symbol_filter)

    def _apply_symbol_filter(self, collection, symbol_filter, create_index_list=False):
        filtered = []
        for idx, item in enumerate(collection):
            if idx not in self._index_to_symbol_map or \
                symbol_filter(self._index_to_symbol_map[idx]):
                filtered.append(idx if create_index_list else item)
        return filtered

    def parse_and_add_profiles(self, profile_properties_list):
        profiles = [self._extract_function(profile_properties)
                    for profile_properties in profile_properties_list]
        mass_values = [self._extract_mass_value(profile_properties)
                       for profile_properties in profile_properties_list]
        self._profiles.append(profiles)
        self._mass_values.append(mass_values)

    def _extract_function(self, profile_properties):
        function_name = ("function=%s," % profile_properties.pop('function'))
        function_props = ["{0}={1}".format(key, value) for key, value in profile_properties.items()]
        function_props = ("%s,%s" % (function_name, (','.join(function_props))))
        return function_props

    def _extract_mass_value(self, profile_properties):
        material_builder = MaterialBuilder()
        mass_value = profile_properties.pop('value', None)
        if mass_value is None:
            symbol = profile_properties.pop('symbol', None)

            if symbol is None:
                raise RuntimeError('Invalid mass specified - ' + str(profile_properties)
                                   + " - either 'value' or 'symbol' must be given.")

            try:
                mass_value = material_builder.setFormula(symbol).build().relativeMolecularMass()
                self._index_to_symbol_map[len(self._mass_values)] = mass_value
            except BaseException as exc:
                raise RuntimeError('Error when parsing mass - ' + str(profile_properties) + ": "
                                   + "\n" + str(exc))
        return mass_value

    def update_from_workspace(self, params_ws, symbol_filter=None):
        function_regex = re.compile("f([0-9]+).([A-z0-9_]+)")
        param_labels = params_ws.getAxis(1).extractValues()
        mass_indices = self._apply_symbol_filter(self._mass_values[0], symbol_filter,
                                                 create_index_list=True)
        num_masses = len(self._mass_values[0])

        for column_idx in range(params_ws.blocksize()):

            for idx, param in enumerate(param_labels):
                if param != "Cost function value":
                    param_re = function_regex.match(param)
                    mass_idx = mass_indices[int(param_re.group(1))]

                    if mass_idx >= num_masses:
                        continue

                    param_name = param_re.group(2).lower()
                    if param_name == 'width' and \
                            isinstance(self._profiles[column_idx][mass_idx].get(param_name, None), list):
                        self._profiles[column_idx][mass_idx][param_name][1] = params_ws.dataY(idx)[column_idx]
                    elif 'symbol' not in self._profiles[column_idx][mass_idx] or param_name != 'mass':
                        self._profiles[column_idx][mass_idx][param_name] = params_ws.dataY(idx)[column_idx]


# -----------------------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------------------

class VesuvioMSHelper(object):

    def __init__(self, beam_radius=2.5, sample_height=5.0, sample_width=5.0, sample_depth=5.0,
                 sample_density=1.0, seed=123456789, num_scatters=3, num_runs=10, num_events=50000,
                 smooth_neighbours=3):
        self._beam_radius = beam_radius
        self._sample_height = sample_height
        self._sample_width = sample_width
        self._sample_depth = sample_depth
        self._sample_density = sample_density
        self._seed = seed
        self._num_scatters = num_scatters
        self._num_runs = num_runs
        self._num_events = num_events
        self._smooth_neighbours = smooth_neighbours
        self._hydrogen_constraints = {}

    def parse_and_add_hydrogen_constraints(self, constraints):
        """
        Parses the specified hydrogen constraints from a list of constraints,
        to a dictionary of the chemical symbols of each constraint mapped to
        their corresponding constraint properties (factor and weight).

        :param constraints: The hydrogen constraints to parse.
        :raise:             A RuntimeError if a constraint doesn't hasn't been
                            given an associated chemical symbol.
        """
        if not isinstance(constraints, dict):

            try:
                for constraint in constraints:
                    self.parse_and_add_hydrogen_constraint(constraint)
            except AttributeError:
                raise RuntimeError("HydrogenConstraints are incorrectly formatted.")

        try:
            self.parse_and_add_hydrogen_constraint(constraints)
        except AttributeError:
            raise RuntimeError("HydrogenConstraints are incorrectly formatted.")

    def parse_and_add_hydrogen_constraint(self, constraint):
        symbol = constraint.pop("symbol", None)

        if symbol is None:
            raise RuntimeError("Invalid hydrogen constraint: " +
                               str(constraint) +
                               " - No symbol provided")
        self._hydrogen_constraints[symbol] = constraint


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

class VesuvioTOFFitRoutine(object):

    def __init__(self, load_helper, ms_helper, mass_profile_collection):
        self._load_helper = load_helper
        self._ms_helper = ms_helper
        self._mass_profile_collection = mass_profile_collection

    def __call__(self, sample_runs, container_runs, spectra,
                 iterations, convergence_threshold):
        if iterations < 1:
            return ValueError('Must perform at least one iteration')
        sample_data = self._load_data(sample_runs)
        container_data = self._load_data(container_runs)

    def _load_data(self, runs, spectra):
        if isinstance(runs, MatrixWorkspace):
            return runs
        else:
            return self._load_helper(runs)

# ----------------------------------------------------------------------------------------

class VesuvioTOFFitRoutineIteration(object):

    def __init__(self, ms_helper, mass_profile_collectiom):
        self._ms_helper = ms_helper
        self._mass_profile_collection = mass_profile_collectiom

    def __call__(self, sample_data, container_data):
