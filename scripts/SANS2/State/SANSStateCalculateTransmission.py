# pylint: disable=too-few-public-methods

"""State describing the calculation of the transmission for SANS reduction."""

import json
from SANS2.State.SANSStateBase import (SANSStateBase, sans_parameters, PositiveIntegerParameter,
                                       BoolParameter, PositiveFloatParameter, ClassTypeParameter,
                                       FloatParameter, DictParameter, StringListParameter, StringParameter,
                                       PositiveFloatWithNoneParameter)
from SANS2.Common.SANSEnumerations import (RebinType, RangeStepType, FitType, DataType,
                                           convert_reduction_data_type_to_string)
from SANS2.Common.SANSConfigurations import SANSConfigurations


# ------------------------------------------------
# SANSStateAdjustment
# ------------------------------------------------
class SANSStateCalculateTransmission(object):
    pass


@sans_parameters
class SANSStateTransmissionFit(SANSStateBase):
    fit_type = ClassTypeParameter(FitType)
    polynomial_order = PositiveIntegerParameter()
    wavelength_low = PositiveFloatWithNoneParameter()
    wavelength_high = PositiveFloatWithNoneParameter()

    def __init__(self):
        super(SANSStateTransmissionFit, self).__init__()
        self.fit_type = FitType.Linear
        self.polynomial_order = 0

    def validate(self):
        is_invalid = {}
        if self.fit_type is not FitType.Polynomial and self.polynomial_order != 0:
            is_invalid.update({"SANSStateTransmissionFit": "Can only set a polynomial order for polynomial"
                                                           " fitting, but selected {0}".format(self.fit_type)})
        if (self.wavelength_low is not None and self.wavelength_high is None) or \
                (self.wavelength_low is None and self.wavelength_high is not None):
            is_invalid.update({"SANSStateTransmissionFit": "Either both an upper and a lower limit "
                                                           "have to be specified or none."})
        if self.wavelength_low is not None and self.wavelength_high is not None:
            if self.wavelength_low > self.wavelength_high:
                is_invalid.update({"SANSStateTransmissionFit": "The lower wavelength limit {0} should be smaller "
                                                               "than the upper wavelength limit {1}"
                                                               ".".format(self.wavelength_low, self.wavelength_high)})
        if is_invalid:
            raise ValueError("SANSStateTransmissionFit: The provided inputs are illegal. "
                             "Please see: {0}".format(json.dumps(is_invalid)))


@sans_parameters
class SANSStateCalculateTransmissionISIS(SANSStateBase, SANSStateCalculateTransmission):
    # -----------------------
    # Transmission
    # -----------------------
    transmission_radius_on_detector = PositiveFloatParameter()
    transmission_roi_files = StringListParameter()
    transmission_mask_files = StringListParameter()

    default_transmission_monitor = PositiveIntegerParameter()
    transmission_monitor = PositiveIntegerParameter()

    default_incident_monitor = PositiveIntegerParameter()
    incident_monitor = PositiveIntegerParameter()

    # ----------------------
    # Prompt peak correction
    # ----------------------
    prompt_peak_correction_min = PositiveFloatParameter()
    prompt_peak_correction_max = PositiveFloatParameter()

    # ----------------
    # Wavelength rebin
    # ----------------
    rebin_type = ClassTypeParameter(RebinType)
    wavelength_low = PositiveFloatParameter()
    wavelength_high = PositiveFloatParameter()
    wavelength_step = PositiveFloatParameter()
    wavelength_step_type = ClassTypeParameter(RangeStepType)

    use_full_wavelength_range = BoolParameter()
    wavelength_full_range_low = PositiveFloatParameter()
    wavelength_full_range_high = PositiveFloatParameter()

    # -----------------------
    # Background correction
    # ----------------------
    background_TOF_general_start = FloatParameter()
    background_TOF_general_stop = FloatParameter()
    background_TOF_monitor_start = DictParameter()
    background_TOF_monitor_stop = DictParameter()
    background_TOF_roi_start = FloatParameter()
    background_TOF_roi_stop = FloatParameter()

    # -----------------------
    # Fit
    # ----------------------
    fit = DictParameter()

    def __init__(self):
        super(SANSStateCalculateTransmissionISIS, self).__init__()
        self.background_TOF_monitor_start = {}
        self.background_TOF_monitor_stop = {}
        self.fit = {convert_reduction_data_type_to_string(DataType.Sample): SANSStateTransmissionFit(),
                    convert_reduction_data_type_to_string(DataType.Can): SANSStateTransmissionFit()}
        self.use_full_wavelength_range = False

    def validate(self):
        is_invalid = {}
        # -----------------
        # Incident monitor
        # -----------------
        if self.incident_monitor is None and self.default_incident_monitor is None:
            is_invalid.update({"incident_monitor": "An incident monitor must be specified."})

        # --------------
        # Transmission, either we need some ROI (ie radius, roi files /mask files) or a transmission monitor
        # --------------
        has_no_transmission_monitor_setting = self.transmission_monitor is None and\
                                              self.default_transmission_monitor is None
        has_no_transmission_roi_setting = self.transmission_radius_on_detector is None and \
                                          self.transmission_roi_files is None
        if has_no_transmission_monitor_setting and has_no_transmission_roi_setting:
            is_invalid.update({"transmission_monitor setting": "A transmission monitor or a transmission radius or"
                                                               " transmission files need to be specified."})

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

        if self.use_full_wavelength_range:
            if self.wavelength_full_range_low is None or self.wavelength_full_range_high is None:
                is_invalid.update({"wavelength_full_range_low": "The full wavelength range is not know."})
            if self.wavelength_full_range_low is not None and self.wavelength_full_range_high is not None and \
                            self.wavelength_full_range_low > self.wavelength_full_range_high:
                is_invalid.update({"wavelength_full_range_low": "The full lower wavelength range is larger than the " \
                                                                "upper wavelength range. The values are {0} and " \
                                                                "{1}".format(self.wavelength_full_range_low,
                                                                             self.wavelength_full_range_high)})
        # ----------------------
        # Background correction
        # ----------------------
        if (self.background_TOF_roi_start is not None and self.background_TOF_roi_stop is None) or \
           (self.background_TOF_roi_start is None and self.background_TOF_roi_stop is not None):
            is_invalid.update({"background_TOF_roi_start": "Only the start or the stop value of the ROI "
                                                           "background TOF correction was specified. Either both "
                                                           "or none has to be specified."})
        if self.background_TOF_roi_start is not None and self.background_TOF_roi_stop is not None:
            if self.background_TOF_roi_start > self.background_TOF_roi_stop:
                is_invalid.update({"background_TOF_roi_start": "The start value of the ROI background TOF "
                                                               "correction is larger than the stop value. The "
                                                               "start value is {0} and the stop value is "
                                                               "{1}".format(self.background_TOF_roi_start,
                                                                            self.background_TOF_roi_stop)})
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
                               "background TOF correction was specified. Either both or none has to be specified."})
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

        # -----
        # Fit
        # -----
        self.fit[convert_reduction_data_type_to_string(DataType.Sample)].validate()
        self.fit[convert_reduction_data_type_to_string(DataType.Can)].validate()

        if is_invalid:
            raise ValueError("SANSStateCalculateTransmissionISIS: The provided inputs are illegal. "
                             "Please see: {0}".format(json.dumps(is_invalid)))


class SANSStateCalculateTransmissionLOQ(SANSStateCalculateTransmissionISIS):
    def __init__(self):
        super(SANSStateCalculateTransmissionLOQ, self).__init__()
        # Set the LOQ full wavelength range
        self.wavelength_full_range_low = SANSConfigurations.LOQ.wavelength_full_range_low
        self.wavelength_full_range_high = SANSConfigurations.LOQ.wavelength_full_range_high

        # Set the LOQ default range for prompt peak correction
        self.prompt_peak_correction_min = SANSConfigurations.LOQ.prompt_peak_correction_min
        self.prompt_peak_correction_max = SANSConfigurations.LOQ.prompt_peak_correction_max

    def validate(self):
        super(SANSStateCalculateTransmissionLOQ, self).validate()


class SANSStateCalculateTransmissionSANS2D(SANSStateCalculateTransmissionISIS):
    def __init__(self):
        super(SANSStateCalculateTransmissionSANS2D, self).__init__()
        # Set the LOQ full wavelength range
        self.wavelength_full_range_low = SANSConfigurations.SANS2D.wavelength_full_range_low
        self.wavelength_full_range_high = SANSConfigurations.SANS2D.wavelength_full_range_high

    def validate(self):
        super(SANSStateCalculateTransmissionSANS2D, self).validate()


class SANSStateCalculateTransmissionLARMOR(SANSStateCalculateTransmissionISIS):
    def __init__(self):
        super(SANSStateCalculateTransmissionLARMOR, self).__init__()
        # Set the LOQ full wavelength range
        self.wavelength_full_range_low = SANSConfigurations.LARMOR.wavelength_full_range_low
        self.wavelength_full_range_high = SANSConfigurations.LARMOR.wavelength_full_range_high

    def validate(self):
        super(SANSStateCalculateTransmissionLARMOR, self).validate()

# -----------------------------------------------
# SANSStateCalculateTransmission setup for other facilities/techniques/scenarios.
# Needs to derive from SANSStateCalculateTransmission and SANSStateBase and fulfill its contract.
# -----------------------------------------------
