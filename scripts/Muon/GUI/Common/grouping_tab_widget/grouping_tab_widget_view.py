from __future__ import (absolute_import, division, print_function)

from PyQt4 import QtCore, QtGui


class GroupingTabView(QtGui.QWidget):

    def __init__(self, grouping_table, pairing_table, parent=None):
        super(GroupingTabView, self).__init__(parent)

        self._grouping_table = grouping_table
        self._pairing_table = pairing_table

        self.setup_interface_layout()

    def setup_interface_layout(self):
        self.setObjectName("GroupingTabView")
        self.resize(500, 500)

        self.vertical_layout = QtGui.QVBoxLayout(self)
        self.vertical_layout.setObjectName("verticalLayout")
        if self._grouping_table:
            self.vertical_layout.addWidget(self._grouping_table)
        if self._pairing_table:
            self.vertical_layout.addWidget(self._pairing_table)

        self.setLayout(self.vertical_layout)

    def set_grouping_table(self, table):
        self._grouping_table = table

    def set_pairing_table(self, table):
        self._pairing_table = table

    def update_tables(self):
        self._grouping_table.update_view_from_model()
        self._pairing_table.update_view_from_model()

    def on_grouping_table_changed(self, slot):
        self._grouping_table.dataChanged.connect(slot)

    def on_pairing_table_changed(self, slot):
        self._pairing_table.dataChanged.connect(slot)