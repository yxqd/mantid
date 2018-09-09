from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict

from Muon.GUI.Common.pairing_table_widget.pairing_table_widget_presenter import MuonPair
from Muon.GUI.Common.grouping_table_widget.grouping_table_widget_presenter import MuonGroup


class GroupingTabModel(object):
    """
    The model for the grouping tab should be shared between all widgets of the tab.
    It keeps a record of the groups and pairs defined for the current instance of the interface.

    pairs and groups should be of type MuonGroup and MuonPair respectively.
    """

    def __init__(self):
        # ordered dictionary to preserve order in which use enters data
        self._groups = OrderedDict()
        self._pairs = OrderedDict()

    @property
    def groups(self):
        return self._groups.values()

    @property
    def pairs(self):
        return self._pairs.values()

    @property
    def group_names(self):
        return self._groups.keys()

    @property
    def pair_names(self):
        return self._pairs.keys()

    @property
    def group_and_pair_names(self):
        return self._groups.keys() + self._pairs.keys()

    def clear_groups(self):
        self._groups = OrderedDict()

    def clear_pairs(self):
        self._pairs = OrderedDict()

    def clear(self):
        self.clear_groups()
        self.clear_pairs()

    def add_group(self, group):
        assert isinstance(group, MuonGroup)
        self._groups[group.name] = group

    def add_pair(self, pair):
        assert isinstance(pair, MuonPair)
        self._pairs[pair.name] = pair

    def remove_groups_by_name(self, name_list):
        for name in name_list:
            del self._groups[name]

    def remove_pairs_by_name(self, name_list):
        for name in name_list:
            del self._pairs[name]

    def construct_empty_group(self, group_index):
        group_index = 0
        new_group_name = "group_" + str(group_index)
        while new_group_name in self.group_names:
            group_index += 1
            new_group_name = "group_" + str(group_index)
        return MuonGroup(group_name=new_group_name, detector_IDs=[1])

    def construct_empty_pair(self, pair_index):
        pair_index = 0
        new_pair_name = "pair_" + str(pair_index)
        while new_pair_name in self.pair_names:
            pair_index += 1
            new_pair_name = "pair_" + str(pair_index)
        group1 = self.group_names[0]
        group2 = self.group_names[1]
        return MuonPair(pair_name=new_pair_name,
                        group1_name=group1, group2_name=group2, alpha=1.0)

    def construct_empty_pair_with_group_names(self, name1, name2):
        print("EMPTY PAIR")
        pair_index = 0
        new_pair_name = "pair_" + str(pair_index)
        while new_pair_name in self.pair_names:
            pair_index += 1
            new_pair_name = "pair_" + str(pair_index)
        return MuonPair(pair_name=new_pair_name,
                        group1_name=name1, group2_name=name2, alpha=1.0)
