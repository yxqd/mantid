from __future__ import (absolute_import, division, print_function)

import six
import re

from Muon.GUI.Common import run_string_utils as run_utils


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

        self._view.on_user_changes_group_name(self.validate_group_name)
        self._view.on_user_changes_detector_IDs(self.validate_detector_IDs)

        self._view.on_table_data_changed(self.handle_data_change)

    def validate_group_name(self, text):
        if sum(text == name for name in self._model.group_names) > 1:
            self._view.warning_popup("Groups must have unique names")
            return False
        if re.match("^\w+$", text):
            return True
        return False

    @staticmethod
    def validate_detector_IDs(text):
        if re.match("^[0-9]*([0-9]+[,-]{0,1})*[0-9]+$", text):
            return True
        return False

    def show(self):
        self._view.show()

    def add_group(self, group):
        self.add_group_to_view(group)
        self._model.add_group(group)
        self._view.notify_data_changed()

    def add_group_to_view(self, group):
        self._view.disable_updates()
        assert isinstance(group, MuonGroup)
        if self._view.num_rows() > 20:
            self._view.warning_popup("Cannot add more than 20 groups.")
            self._view.enable_updates()
            return

        #entry = [str(group.name), detector_list_to_string(group.detectors), str(group.n_detectors)]
        entry = [str(group.name), run_utils.run_list_to_string(group.detectors), str(group.n_detectors)]
        self._view.add_entry_to_table(entry)
        self._view.enable_updates()

    def handle_add_group_button_clicked(self):
        group = self._model.construct_empty_group(self._view.num_rows() + 1)
        self.add_group(group)

    def handle_remove_group_button_clicked(self):
        group_names = self._view.get_selected_group_names()
        if not group_names:
            self._view.remove_last_row()
        else:
            self._view.remove_selected_groups()
            self._model.remove_groups_by_name(group_names)
        self._view.notify_data_changed()

    def handle_data_change(self):
        self.update_model_from_view()
        self.update_view_from_model()
        self._view.notify_data_changed()

    def update_model_from_view(self):
        table = self._view.get_table_contents()
        print("Model update from view : ", table)
        self._model.clear_groups()
        for entry in table:
            # TODO parse string to list of detectors
            detector_list = run_utils.run_string_to_list(str(entry[1]))
            group = MuonGroup(group_name=str(entry[0]), detector_IDs=detector_list)
            self._model.add_group(group)

    def update_view_from_model(self):
        self._view.disable_updates()

        self._view.clear()
        for group in self._model.groups.values():
            self.add_group_to_view(group)

        self._view.enable_updates()
