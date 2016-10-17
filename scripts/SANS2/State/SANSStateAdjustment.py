# pylint: disable=too-few-public-methods

"""State describing the adjustment workspace creation of the SANS reduction."""

import json
from SANS2.State.StateFunctions import validator_sub_state
from SANS2.State.SANSStateBase import (SANSStateBase, TypedParameter, sans_parameters)
from SANS2.State.SANSStateCalcualteTransmission import SANSStateCalcualteTransmission
from SANS2.State.SANSStateNormalizeToMonitor import SANSStateNormalizeToMonitor
from SANS2.State.SANSStateWavelengthAndPixelAdjustment import SANSStateWavelengthAndPixelAdjustment


# ------------------------------------------------
# SANSStateAdjustment
# ------------------------------------------------
class SANSStateAdjustment(object):
    pass


@sans_parameters
class SANSStateAdjustmentISIS(SANSStateBase, SANSStateAdjustment):
    calculate_transmission = TypedParameter(SANSStateCalcualteTransmission, validator_sub_state)
    normalize_to_monitor = TypedParameter(SANSStateNormalizeToMonitor, validator_sub_state)
    wavelength_and_pixel_adjustment = TypedParameter(SANSStateWavelengthAndPixelAdjustment, validator_sub_state)

    def __init__(self):
        super(SANSStateAdjustmentISIS, self).__init__()

    def validate(self):
        is_invalid = {}


        if is_invalid:
            raise ValueError("SANSStateAdjustmentISIS: The provided inputs are illegal. "
                             "Please see: {0}".format(json.dumps(is_invalid)))

# -----------------------------------------------
# SANSStateAdjustment setup for other facilities/techniques/scenarios.
# Needs to derive from SANSStateAdjustment and SANSStateBase and fulfill its contract.
# -----------------------------------------------
