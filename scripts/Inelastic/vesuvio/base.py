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

class VesuvioLoadHelper(object):
    """
    A helper class for loading Vesuvio data from the input of a user script.
    """

    def __init__(self, diff_mode, fit_mode, param_file):
        self._diff_mode = diff_mode
        self._fit_mode = fit_mode
        self._param_file = param_file
        self._instrument = VESUVIO()

    def load_and_crop_runs(self, runs, spectra, rebin_params=None):
        sum_spectra = (self._fit_mode == 'bank')
        loaded = self.load_runs(runs, spectra, sum_spectra)
        cropped = self._crop_workspace(loaded)

        if rebin_params is not None:
            return self._rebin_workspace(cropped, rebin_params)
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

    def _rebin_workspace(self, workspace, rebin_params):
        rebin_alg = AlgorithmManager.createUnmanaged("Rebin")
        rebin_alg.setProperty("InputWorkspace", workspace)
        rebin_alg.setProperty("Params", rebin_params)
        rebin_alg.setProperty("OutputWorkspace", "rebinned")
        rebin_alg.execute()
        return rebin_alg.getProperty("OutputWorkspace")

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
