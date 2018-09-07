from __future__ import (absolute_import, division, print_function)

import six


def detector_list_to_string(detector_list):
    return ",".join([str(i) for i in detector_list])


class MuonGroup:
    """Simple struct to store information on a detector group.

    The name is set at initialization and after that cannot be changed.
    The detector list can be modified by passing a list of ints (type checks for this)
    The number of detetors is stored
    """

    def __init__(self, group_name="", detector_IDs=[]):
        self._group_name = group_name

        self._detector_IDs = None
        self.detectors = detector_IDs

    @property
    def name(self):
        return self._group_name

    @property
    def detectors(self):
        return self._detector_IDs

    @property
    def n_detectors(self):
        return len(self.detectors)

    @detectors.setter
    def detectors(self, detector_IDs):
        if isinstance(detector_IDs, six.string_types):
            raise ValueError("detectors must be a list of ints.")
        elif isinstance(detector_IDs, list):
            if sum([not isinstance(item, int) for item in detector_IDs]) == 0:
                self._detector_IDs = detector_IDs
            else:
                raise ValueError("detectors must be a list of ints.")
        else:
            raise ValueError("detectors must be a list of ints.")


class GroupingTablePresenter(object):

    def __init__(self, view, model):
        self._view = view
        self._model = model

        self._view.on_add_group_button_clicked(self.handle_add_group_button_clicked)
        self._view.on_remove_group_button_clicked(self.handle_remove_group_button_clicked)

    def show(self):
        self._view.show()

    def add_group(self, group):
        assert isinstance(group, MuonGroup)
        if self._view.num_rows() > 20:
            self._view.warning_popup("Cannot add more than 20 groups.")
            return

        entry = [str(group.name), detector_list_to_string(group.detectors), str(group.n_detectors)]
        self._view.add_entry_to_table(entry)
        self._model.add_group(group)

    def handle_add_group_button_clicked(self):
        group = self._model.construct_empty_group(self._view.num_rows() + 1)
        self.add_group(group)

    def handle_remove_group_button_clicked(self):
        group_names = self._view.get_selected_group_names()
        print(group_names)
        if not group_names:
            self._view.remove_last_row()
        else:
            self._view.remove_selected_groups()
            self._model.remove_groups_by_name(group_names)
