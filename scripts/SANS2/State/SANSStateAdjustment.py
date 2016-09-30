# pylint: disable=too-few-public-methods

"""State describing the adjustment workspace creation of the SANS reduction."""

import json
from SANS2.State.SANSStateBase import (SANSStateBase, sans_parameters, PositiveIntegerParameter,
                                       BoolParameter, PositiveFloatParameter, ClassTypeParameter,
                                       FloatParameter, DictParameter, StringListParameter, StringParmeter)
from SANS2.Common.SANSEnumerations import (RebinType, RangeStepType)
)


# ------------------------------------------------
# SANSStateAdjustment
# ------------------------------------------------
class SANSStateAdjustment(object):
    pass


@sans_parameters
class SANSStateAdjustmentISIS(SANSStateBase, SANSStateAdjustment):
    # -------------------------------------
    # General parameters
    # -------------------------------------
    use_prompt_peak_correction = BoolParameter()
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

    # -------------------------------------
    # Parameters for Monitor Normalization
    # -------------------------------------
    incident_monitor = PositiveIntegerParameter()

    # ---------------------------------------
    # Parameters for Transmission Calculation
    # ---------------------------------------
    use_full_wavelength_range = BoolParameter()
    wavelength_range_low = PositiveFloatParameter()
    wavelength_range_high = PositiveFloatParameter()
    default_wavelength_range_low = PositiveFloatParameter()
    default_wavelength_range_high = PositiveFloatParameter()

    background_TOF_ROI_start = FloatParameter()
    background_TOF_ROI_stop = FloatParameter()

    transmission_radius_on_detector = PositiveFloatParameter()
    transmission_roi_files = StringListParameter()
    transmission_mask_files = StringListParameter()
    incident_monitor_for_transmission_calculation = PositiveIntegerParameter()

    default_transmission_monitor = PositiveIntegerParameter()
    transmission_monitor = PositiveIntegerParameter()

    # ---------------------------------------
    # Parameters for adjustments
    # ---------------------------------------
    pixel_adjustment_file = StringParmeter()
    wavelength_adjustment_file = StringParmeter()

    def __init__(self):
        super(SANSStateAdjustmentISIS, self).__init__()

    def validate(self):
        is_invalid = {}
        if is_invalid:
            raise ValueError("SANSStateMoveDetectorISIS: The provided inputs are illegal. "
                             "Please see: {0}".format(json.dumps(is_invalid)))

# -----------------------------------------------
# SANSStateAdjustment setup for other facilities/techniques/scenarios.
# Needs to derive from SANSStateAdjustment and SANSStateBase and fulfill its contract.
# -----------------------------------------------
