from __future__ import (absolute_import, division, print_function)


class GroupingTabPresenter(object):

    def __init__(self, view, model,
                 grouping_table_widget=None,
                 pairing_table_widget=None):
        self._view = view
        self._model = model

        self.grouping_table_widget = grouping_table_widget
        self.pairing_table_widget = pairing_table_widget

        self._view.on_grouping_table_changed(self.pairing_table_widget.update_view_from_model)
        self._view.on_pairing_table_changed(self.grouping_table_widget.update_view_from_model)

    def show(self):
        self._view.show()

    def update_tables(self):
        self.grouping_table_widget.update_view_from_model()
        self.pairing_table_widget.update_view_from_model()
