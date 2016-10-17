# pylint: disable=too-few-public-methods

"""State describing the creation of pixel and wavelength adjustment workspaces for SANS reduction."""

import json
from SANS2.State.SANSStateBase import (SANSStateBase, sans_parameters, StringParameter,
                                       ClassTypeParameter, PositiveFloatParameter)
from SANS2.Common.SANSEnumerations import (RangeStepType)


# ------------------------------------------------
# SANSStateAdjustment
# ------------------------------------------------
class SANSStateWavelengthAndPixelAdjustment(object):
    pass


@sans_parameters
class SANSStateWavelengthAndPixelAdjustmentISIS(SANSStateBase, SANSStateWavelengthAndPixelAdjustment):
    pixel_adjustment_file = StringParameter()
    wavelength_adjustment_file = StringParameter()

    wavelength_low = PositiveFloatParameter()
    wavelength_high = PositiveFloatParameter()
    wavelength_step = PositiveFloatParameter()
    wavelength_step_type = ClassTypeParameter(RangeStepType)

    def __init__(self):
        super(SANSStateWavelengthAndPixelAdjustmentISIS, self).__init__()

    def validate(self):
        is_invalid = {}

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
        if is_invalid:
            raise ValueError("SANSStateWavelengthAndPixelAdjustmentISIS: The provided inputs are illegal. "
                             "Please see: {0}".format(json.dumps(is_invalid)))


# -----------------------------------------------
# SANSStateNormalizeMonitor setup for other facilities/techniques/scenarios.
# Needs to derive from SANSStateNormalizeMonitor and SANSStateBase and fulfill its contract.
# -----------------------------------------------
