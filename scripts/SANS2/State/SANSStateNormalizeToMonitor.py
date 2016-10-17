# pylint: disable=too-few-public-methods

"""State describing the normalization to the incident monitor for SANS reduction."""

import json
from SANS2.State.SANSStateBase import (SANSStateBase, sans_parameters, PositiveIntegerParameter,
                                       PositiveFloatParameter, FloatParameter, ClassTypeParameter, DictParameter)
from SANS2.Common.SANSEnumerations import (RebinType, RangeStepType)


# ------------------------------------------------
# SANSStateAdjustment
# ------------------------------------------------
class SANSStateNormalizeToMonitor(object):
    pass


@sans_parameters
class SANSStateNormalizeToMonitorISIS(SANSStateBase, SANSStateNormalizeToMonitor):
    prompt_peak_correction_min = PositiveFloatParameter()
    prompt_peak_correction_max = PositiveFloatParameter()

    rebin_type = ClassTypeParameter(RebinType)
    wavelength_low = PositiveFloatParameter()
    wavelength_high = PositiveFloatParameter()
    wavelength_step = PositiveFloatParameter()
    wavelength_step_type = ClassTypeParameter(RangeStepType)

    background_TOF_general_start = FloatParameter()
    background_TOF_general_stop = FloatParameter()
    background_TOF_monitor_start = DictParameter()
    background_TOF_monitor_stop = DictParameter()

    incident_monitor = PositiveIntegerParameter()

    def __init__(self):
        super(SANSStateNormalizeToMonitorISIS, self).__init__()
        self.background_TOF_monitor_start = {}
        self.background_TOF_monitor_stop = {}

    def validate(self):
        is_invalid = {}
        # -----------------
        # incident Monitor
        # -----------------
        if self.incident_monitor is None:
            is_invalid.update({"incident_monitor": "An incident monitor must be specified."})

        # -----------------
        # Prompt peak
        # -----------------
        if (self.prompt_peak_correction_min is None and self.prompt_peak_correction_max is not None) or \
                (self.prompt_peak_correction_min is not None and self.prompt_peak_correction_max is None):
            is_invalid.update({"prompt_peak_correction_min": "You have to specify either the start and stop value for"
                                                             " the prompt peak correction or none."})
        if self.prompt_peak_correction_min is not None and self.prompt_peak_correction_max is not None:
            if self.prompt_peak_correction_min > self.prompt_peak_correction_max:
                is_invalid.update({"prompt_peak_correction_min": "The start value for the prompt peak correction needs "
                                                                 "to be smaller than the stop value. "
                                                                 "The start value was {0} and the stop value "
                                                                 "{1}.".format(self.prompt_peak_correction_min,
                                                                               self.prompt_peak_correction_max)})

        # -----------------
        # Wavelength rebin
        # -----------------
        if self.rebin_type is None:
            is_invalid.update({"rebin_type": "A rebin type has to be specified."})
        if self.wavelength_step_type is None:
            is_invalid.update({"wavelength_step_type": "A wavelength range step type has to be specified."})
        if self.wavelength_step is None:
            is_invalid.update({"wavelength_step": "A wavelength step has to be specified."})
        if self.wavelength_low is None:
            is_invalid.update({"wavelength_low": "A lower wavelength value for rebinning has to be specified."})
        if self.wavelength_high is None:
            is_invalid.update({"wavelength_high": "An high wavelength value for rebinning has to be specified."})
        if self.wavelength_low is not None and self.wavelength_high is not None:
            if self.wavelength_low > self.wavelength_high:
                is_invalid.update({"wavelength_high": "The lower wavelength bound needs to be smaller than the upper "
                                                      "bound. The lower bound is {0} and the upper "
                                                      "is {1}.".format(self.wavelength_low, self.wavelength_high)})

        # ----------------------
        # Background correction
        # ----------------------
        if (self.background_TOF_general_start is not None and self.background_TOF_general_stop is None) or \
           (self.background_TOF_general_start is None and self.background_TOF_general_stop is not None):
            is_invalid.update({"background_TOF_general_start": "Only the start or the stop value of the general "
                                                               "background TOF correction was specified. Either both "
                                                               "or none has to be specified."})
        if self.background_TOF_general_start is not None and self.background_TOF_general_stop is not None:
            if self.background_TOF_general_start > self.background_TOF_general_stop:
                is_invalid.update({"background_TOF_general_start": "The start value of the general background TOF "
                                                                   "correction is larger than the stop value. The "
                                                                   "start value is {0} and the stop value is "
                                                                   "{1}".format(self.background_TOF_general_start,
                                                                                self.background_TOF_general_stop)})
        if (self.background_TOF_monitor_start is not None and self.background_TOF_monitor_stop is None) or \
           (self.background_TOF_monitor_start is None and self.background_TOF_monitor_stop is not None):
            is_invalid.update({"background_TOF_monitor_start value": "Only the start or the stop value of the monitor "
                                                               "background TOF correction was specified. Either both "
                                                               "or none has to be specified."})
        if self.background_TOF_monitor_start is not None and self.background_TOF_monitor_stop is not None:
            if len(self.background_TOF_monitor_start) != len(self.background_TOF_monitor_stop):
                is_invalid.update({"background_TOF_monitor_start length": "The number of entries between for start and"
                                                                          " stop values of monitor backgrounds don't "
                                                                          "match. There are {0} start values and {1} "
                "stop values".format(len(self.background_TOF_monitor_start), len(self.background_TOF_monitor_stop))})
            for key_start, value_start in self.background_TOF_monitor_start.items():
                if key_start not in self.background_TOF_monitor_stop:
                    is_invalid.update({"background_TOF_monitor_start key": "Monitor {0} could not be found for both the"
                                                                           " start and stop value.".format(key_start)})
                else:
                    value_stop = self.background_TOF_monitor_stop[key_start]
                    if value_start > value_stop:
                        is_invalid.update({"background_TOF_monitor_start value": "The start value {0} is larger than "
                                                                                 "the stop value {1} for monitor "
                                                                                 "{2}.".format(value_start, value_stop,
                                                                                               key_start)})

        if is_invalid:
            raise ValueError("SANSStateMoveDetectorISIS: The provided inputs are illegal. "
                             "Please see: {0}".format(json.dumps(is_invalid)))


@sans_parameters
class SANSStateNormalizeToMonitorLOQ(SANSStateNormalizeToMonitorISIS):
    def __init__(self):
        super(SANSStateNormalizeToMonitorLOQ, self).__init__()
        # Set the LOQ default range for prompt peak correction
        self.prompt_peak_correction_min = 19000.0
        self.prompt_peak_correction_max = 20500.0

    def validate(self):
        super(SANSStateNormalizeToMonitorLOQ, self).validate()

# -----------------------------------------------
# SANSStateNormalizeMonitor setup for other facilities/techniques/scenarios.
# Needs to derive from SANSStateNormalizeMonitor and SANSStateBase and fulfill its contract.
# -----------------------------------------------
