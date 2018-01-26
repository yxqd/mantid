from sans.common.constants import EMPTY_NAME
from sans.common.general_functions import create_unmanaged_algorithm

def perform_load(states):
    for state in states:
        serialized_state = state.property_manager
        load_algorithm = create_load_algorithm(serialized_state)
        load_algorithm.execute()

def create_load_algorithm(serialized_state):
    load_name = "SANSLoad"
    load_options = {"SANSState": serialized_state,
                    "PublishToCache": True,
                    "UseCached": True,
                    "SampleScatterWorkspace": EMPTY_NAME,
                    "SampleScatterMonitorWorkspace": EMPTY_NAME,
                    "SampleTransmissionWorkspace": EMPTY_NAME,
                    "SampleDirectWorkspace": EMPTY_NAME,
                    "CanScatterWorkspace": EMPTY_NAME,
                    "CanScatterMonitorWorkspace": EMPTY_NAME,
                    "CanTransmissionWorkspace": EMPTY_NAME,
                    "CanDirectWorkspace": EMPTY_NAME}
    return create_unmanaged_algorithm(load_name, **load_options)