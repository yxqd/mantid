from __future__ import (absolute_import, division, print_function)

from Muon.GUI.Common.grouping_table_widget.grouping_table_widget_presenter import MuonGroup


class GroupingTableModel(object):

    def __init__(self):
        self._groups = {}

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

    def clear(self):
        self._groups = {}
