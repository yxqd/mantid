from __future__ import (absolute_import, division, print_function)


class GroupingTabPresenter(object):
    """

    The grouping tab presenter is responsible for synchronizing the group and pair tables.
    """

    def __init__(self, view, model,
                 grouping_table_widget=None,
                 pairing_table_widget=None):
        self._view = view
        self._model = model

        self.grouping_table_widget = grouping_table_widget
        self.pairing_table_widget = pairing_table_widget

        self._view.on_clear_grouping_button_clicked(self.on_clear_requested)

        # Synchronize the two tables
        self._view.on_grouping_table_changed(self.pairing_table_widget.update_view_from_model)
        self._view.on_pairing_table_changed(self.grouping_table_widget.update_view_from_model)

        self._view.on_add_pair_requested(self.add_pair_from_grouping_table)

        self._view.set_description_text(self.text_for_description())


    def add_pair_from_grouping_table(self, name1, name2):
        """If user requests to add a pair from the grouping table."""
        pair = self._model.construct_empty_pair_with_group_names(name1, name2)
        self._model.add_pair(pair)
        self.pairing_table_widget.update_view_from_model()

    def text_for_description(self):
        # TODO :  implement this
        text = "EMU longitudinal (?? detectors)"
        return text

    def show(self):
        self._view.show()

    def on_clear_requested(self):
        self._model.clear()
        self.grouping_table_widget.update_view_from_model()
        self.pairing_table_widget.update_view_from_model()