from __future__ import (absolute_import, division, print_function)

import six
import re


def detector_list_to_string(detector_list):
    return ",".join([str(i) for i in detector_list])


class MuonPair:
    """Simple struct to store information on a detector group pair.

    The name is set at initialization and after that cannot be changed.
    The detector list can be modified by passing a list of ints (type checks for this)
    The number of detetors is stored
    """

    def __init__(self, pair_name="", group1_name="", group2_name="", alpha=1.0):
        self._pair_name = pair_name
        self._group1_name = group1_name
        self._group2_name = group2_name
        self._alpha = alpha

    @property
    def name(self):
        return self._pair_name

    @property
    def group1(self):
        return self._group1_name

    @property
    def group2(self):
        return self._group2_name

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, new_alpha):
        if new_alpha >= 0.0:
            self._alpha = new_alpha
        else:
            raise ValueError("Alpha must be > 0.0.")


class PairingTablePresenter(object):

    def __init__(self, view, model):
        self._view = view
        self._model = model

        self._view.on_add_pair_button_clicked(self.handle_add_pair_button_clicked)
        self._view.on_remove_pair_button_clicked(self.handle_remove_pair_button_clicked)

        self._view.on_user_changes_pair_name(self.validate_pair_name)

        self._view.on_table_data_changed(self.handle_data_change)

    def validate_pair_name(self, text):
        if sum(text == name for name in self._model.group_and_pair_names) > 0:
            self._view.warning_popup("Groups and pairs must have unique names")
            return False
        if not re.match("^\w+$", text):
            self._view.warning_popup("Group names should only contain digits, characters and _")
            return False
        return True

    @staticmethod
    def validate_detector_IDs(text):
        if re.match("^[0-9]*([0-9]+[,-]{0,1})*[0-9]+$", text):
            return True
        return False

    def show(self):
        self._view.show()

    def add_pair(self, pair):
        self.add_pair_to_view(pair)
        self._model.add_pair(pair)

    def add_pair_to_view(self, pair):
        self._view.disable_updates()
        self.update_group_selections()
        assert isinstance(pair, MuonPair)
        if self._view.num_rows() > 20:
            self._view.warning_popup("Cannot add more than 20 pairs.")
            self._view.enable_updates()
            return
        entry = [str(pair.name), str(pair.group1), str(pair.group2), str(pair.alpha)]
        self._view.add_entry_to_table(entry)
        self._view.enable_updates()

    def handle_add_pair_button_clicked(self):
        pair = self._model.construct_empty_pair(self._view.num_rows() + 1)
        self.add_pair(pair)

    def handle_remove_pair_button_clicked(self):
        pair_names = self._view.get_selected_pair_names()
        if not pair_names:
            self.remove_last_row_in_view_and_model()
        else:
            self._view.remove_selected_pairs()
            self._model.remove_pairs_by_name(pair_names)

    def remove_last_row_in_view_and_model(self):
        name = self._view.get_table_contents()[-1][0]
        self._view.remove_last_row()
        self._model.remove_pairs_by_name([name])

    def handle_data_change(self):
        self.update_model_from_view()
        self.update_view_from_model()

    def update_model_from_view(self):
        table = self._view.get_table_contents()
        self._model.clear_pairs()
        for entry in table:
            pair = MuonPair(pair_name=str(entry[0]), group1_name=str(entry[1]), group2_name=str(entry[2]),
                            alpha=float(entry[3]))
            self._model.add_pair(pair)

    def update_view_from_model(self):
        self._view.disable_updates()

        self._view.clear()
        self.update_group_selections()
        print(self._view._group_selections)
        for pair in self._model.pairs:
            print(pair)
            self.add_pair_to_view(pair)

        self._view.enable_updates()

    def update_group_selections(self):
        groups = self._model.group_names
        self._view.update_group_selections(groups)
