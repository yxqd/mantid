from __future__ import (absolute_import, division, print_function)

from Muon.GUI.Common.muon_pair import MuonPair


class GroupingTableModel(object):

    def __init__(self):
        self._pairs = {}

    @property
    def pairs(self):
        return self._pairs

    @property
    def pair_names(self):
        return self._pairs.keys()

    def add_pair(self, pair):
        self._pairs[pair.name] = pair

    def construct_empty_pair(self, pair_index):
        return MuonPair(pair_name="pair_" + str(pair_index))

    def remove_pairs_by_name(self, name_list):
        for name in name_list:
            del self._pairs[name]

    def clear(self):
        self._pairs = {}
