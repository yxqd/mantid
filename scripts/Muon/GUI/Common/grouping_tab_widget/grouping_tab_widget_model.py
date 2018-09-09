from __future__ import (absolute_import, division, print_function)

from Muon.GUI.Common.pairing_table_widget.pairing_table_widget_presenter import MuonPair
from Muon.GUI.Common.grouping_table_widget.grouping_table_widget_presenter import MuonGroup

class GroupingTabModel(object):
    """Shared between all widgets of the tab"""

    def __init__(self):
        self._groups = {}
        self._pairs = {}

    @property
    def pairs(self):
        return self._pairs

    @property
    def pair_names(self):
        return self._pairs.keys()

    def add_pair(self, pair):
        print("Add pair : ", pair)
        self._pairs[pair.name] = pair

    def construct_empty_pair(self, pair_index):
        return MuonPair(pair_name="pair_" + str(pair_index))

    def remove_pairs_by_name(self, name_list):
        for name in name_list:
            del self._pairs[name]

    def clear_groups(self):
        self._groups = {}

    def clear_pairs(self):
        self._pairs = {}

    def clear(self):
        self.clear_groups()
        self.clear_pairs()

    @property
    def groups(self):
        return self._groups

    @property
    def group_names(self):
        return self._groups.keys()

    def add_group(self, group):
        self._groups[group.name] = group

    def construct_empty_group(self, group_index):
        return MuonGroup(group_name="group_" + str(group_index), detector_IDs=[1])

    def remove_groups_by_name(self, name_list):
        for name in name_list:
            del self._groups[name]


