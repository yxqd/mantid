# pylint: disable=too-few-public-methods

import copy

from SANS2.State.SANSStateCalculateTransmission import (SANSStateCalculateTransmissionISIS)
from SANS2.Common.SANSEnumerations import SANSInstrument
from SANS2.State.StateBuilder.AutomaticSetters import (automatic_setters)
from SANS2.Common.XMLParsing import get_named_elements_from_ipf_file
from SANS2.Common.SANSFileInformation import (get_instrument_paths_for_sans_file)


# -------------------------------------
# Free functions
# -------------------------------------
def set_default_monitors(normalize_monitor_info, data_info):
    """
    The default incident monitor is stored on the IPF.
    :param normalize_monitor_info: a SANSStateNormalizeMonitor object on which we set the default value
    :param data_info: a SANSStateData object
    """
    file_name = data_info.sample_scatter
    _, ipf_path = get_instrument_paths_for_sans_file(file_name)
    incident_tag = "INCIDENT"
    transmission_tag = "TRANSMISSION"
    monitors_to_search = {incident_tag: "default-transmission-monitor-spectrum",
                          transmission_tag: "default-transmission-monitor-spectrum"}
    found_monitor_spectrum = get_named_elements_from_ipf_file(ipf_path, monitors_to_search, int)
    if incident_tag in found_monitor_spectrum:
        normalize_monitor_info.default_incident_monitor = found_monitor_spectrum[incident_tag]
    if transmission_tag in found_monitor_spectrum:
        normalize_monitor_info.default_transmission_monitor = found_monitor_spectrum[transmission_tag]


# ---------------------------------------
# State builders
# ---------------------------------------
class SANSStateCalculateTransmissionBuilderISIS(object):
    @automatic_setters(SANSStateCalculateTransmissionISIS)
    def __init__(self, data_info):
        super(SANSStateCalculateTransmissionBuilderISIS, self).__init__()
        self._data = data_info
        self.state = SANSStateCalculateTransmissionISIS()
        set_default_monitors(self.state, self._data)

    def build(self):
        self.state.validate()
        return copy.copy(self.state)


class SANSStateCalculateTransmissionBuilderLARMOR(SANSStateCalculateTransmissionBuilderISIS):
    @automatic_setters(SANSStateCalculateTransmissionISIS)
    def __init__(self, data_info):
        super(SANSStateCalculateTransmissionBuilderISIS, self).__init__()
        self._data = data_info
        self.state = SANSStateCalculateTransmissionISIS()
        set_default_monitors(self.state, self._data)

    def build(self):
        self.state.validate()
        return copy.copy(self.state)


# ------------------------------------------
# Factory method for SANSStateNormalizeToMonitorBuilder
# ------------------------------------------
def get_normalize_to_monitor_builder(data_info):
    instrument = data_info.instrument
    if instrument is SANSInstrument.LARMOR or instrument is SANSInstrument.LOQ or instrument is SANSInstrument.SANS2D:
        return SANSStateCalculateTransmissionBuilderISIS(data_info)
    else:
        raise NotImplementedError("SANSStateNormalizeToMonitorBuilder: Could not find any valid normalize to monitor "
                                  "builder for the specified SANSStateData object {0}".format(str(data_info)))
