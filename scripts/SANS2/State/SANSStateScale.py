""" Defines the state of the goemetry and unit scaling."""

import json
from SANS2.State.SANSStateBase import (SANSStateBase, sans_parameters, PositiveFloatParameter, ClassTypeParameter)
from SANS2.Common.SANSEnumerations import SampleShape


# ------------------------------------------------
# SANSStateScale
# ------------------------------------------------
class SANSStateScale(object):
    pass


# -----------------------------------------------
#  SANSStateScale for ISIS
# -----------------------------------------------
@sans_parameters
class SANSStateScaleISIS(SANSStateScale, SANSStateBase):
    shape = ClassTypeParameter(SampleShape)
    thickness = PositiveFloatParameter()
    width = PositiveFloatParameter()
    height = PositiveFloatParameter()
    scale = PositiveFloatParameter()

    def __init__(self):
        super(SANSStateScaleISIS, self).__init__()

    def validate(self):
        is_invalid = dict()

        if is_invalid:
            raise ValueError("SANSStateScale: The provided inputs are illegal. "
                             "Please see: {}".format(json.dumps(is_invalid)))



# -----------------------------------------------
# SANSStateScale setup for other facilities/techniques/scenarios.
# Needs to derive from SANSStateScale and SANSStateBase and fulfill its contract.
# -----------------------------------------------
