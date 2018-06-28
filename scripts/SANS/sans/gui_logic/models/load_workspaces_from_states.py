from __future__ import (absolute_import, division, print_function)

import copy

from sans.common.enums import (SANSDataType)
from sans.algorithm_detail.batch_execution import provide_loaded_data


def load_workspaces(states):
    workspace_to_name = {SANSDataType.SampleScatter: "SampleScatterWorkspace",
                         SANSDataType.SampleTransmission: "SampleTransmissionWorkspace",
                         SANSDataType.SampleDirect: "SampleDirectWorkspace",
                         SANSDataType.CanScatter: "CanScatterWorkspace",
                         SANSDataType.CanTransmission: "CanTransmissionWorkspace",
                         SANSDataType.CanDirect: "CanDirectWorkspace"}

    workspace_to_monitor = {SANSDataType.SampleScatter: "SampleScatterMonitorWorkspace",
                            SANSDataType.CanScatter: "CanScatterMonitorWorkspace"}

    states_copy = copy.copy(states)

    for key, state in states_copy.items():
        workspaces, monitors = provide_loaded_data(state, True, workspace_to_name, workspace_to_monitor)
    return True
