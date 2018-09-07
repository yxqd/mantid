from __future__ import (absolute_import, division, print_function)

from PyQt4 import QtCore, QtGui

from Muon.GUI.Common import table_utils
from Muon.GUI.Common import message_box


class GroupingTableView(QtGui.QWidget):

    def __init__(self, parent=None):
        super(GroupingTableView, self).__init__(parent)

        self.grouping_table = QtGui.QTableWidget(self)
        self.set_up_table()

        self.setup_interface_layout()

    def setup_interface_layout(self):
        self.setObjectName("GroupingTableView")
        self.resize(500, 500)

        self.add_group_button = QtGui.QToolButton()
        self.remove_group_button = QtGui.QToolButton()

        size_policy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        size_policy.setHorizontalStretch(0)
        size_policy.setVerticalStretch(0)
        size_policy.setHeightForWidth(self.add_group_button.sizePolicy().hasHeightForWidth())
        size_policy.setHeightForWidth(self.remove_group_button.sizePolicy().hasHeightForWidth())

        self.add_group_button.setSizePolicy(size_policy)
        self.add_group_button.setMinimumSize(QtCore.QSize(25, 25))
        self.add_group_button.setObjectName("addGroupButton")
        self.add_group_button.setText("+")

        self.remove_group_button.setSizePolicy(size_policy)
        self.remove_group_button.setMinimumSize(QtCore.QSize(25, 25))
        self.remove_group_button.setObjectName("removeGroupButton")
        self.remove_group_button.setText("-")

        self.vertical_layout = QtGui.QVBoxLayout(self)
        self.vertical_layout.setObjectName("horizontalLayout")
        self.vertical_layout.addWidget(self.grouping_table)
        self.vertical_layout.addWidget(self.add_group_button)
        self.vertical_layout.addWidget(self.remove_group_button)
        self.setLayout(self.vertical_layout)

    def set_up_table(self):
        self.grouping_table.setColumnCount(3)
        self.grouping_table.setHorizontalHeaderLabels(QtCore.QString("Group Name;Detector IDs;N Detectors").split(";"))
        header = self.grouping_table.horizontalHeader()
        header.setResizeMode(0, QtGui.QHeaderView.Stretch)
        header.setResizeMode(1, QtGui.QHeaderView.Stretch)
        header.setResizeMode(2, QtGui.QHeaderView.ResizeToContents)
        # table_utils.setTableHeaders(self.grouping_table)
        vertical_headers = self.grouping_table.verticalHeader()
        vertical_headers.setMovable(True)
        vertical_headers.setResizeMode(QtGui.QHeaderView.ResizeToContents)

    def _detector_cell_widget(self):
        """
        Create a regex protected LineEdit for the detector entry
        """
        edit = QtGui.QLineEdit(self)
        run_string_regex = "^[0-9]*([0-9]+[,-]{0,1})*[0-9]+$"
        regex = QtCore.QRegExp(run_string_regex)
        validator = QtGui.QRegExpValidator(regex)
        edit.setValidator(validator)
        edit.setStyleSheet("border: 0px solid white;")
        return edit

    def _group_name_cell_widget(self):
        edit = QtGui.QLineEdit(self)
        run_string_regex = "[0-9,a-z,A-Z]*"
        regex = QtCore.QRegExp(run_string_regex)
        validator = QtGui.QRegExpValidator(regex)
        edit.setValidator(validator)
        edit.setStyleSheet("border: 0px solid white;")
        return edit

    def add_entry_to_table(self, row_entries):
        assert len(row_entries) == self.grouping_table.columnCount()

        row_position = self.grouping_table.rowCount()
        self.grouping_table.insertRow(row_position)
        for i, entry in enumerate(row_entries):
            item = QtGui.QTableWidgetItem(entry)
            if i == 0:
                group_name_widget = self._group_name_cell_widget()
                group_name_widget.setText(entry)
                self.grouping_table.setCellWidget(row_position, 0, group_name_widget)
                continue
            if i == 1:
                detector_widget = self._detector_cell_widget()
                detector_widget.setText(entry)
                self.grouping_table.setCellWidget(row_position, 1, detector_widget)
                continue
            if i == 2:
                item.setFlags(QtCore.Qt.ItemIsEnabled)
            self.grouping_table.setItem(row_position, i, item)

    def on_add_group_button_clicked(self, slot):
        self.add_group_button.clicked.connect(slot)

    def on_remove_group_button_clicked(self, slot):
        self.remove_group_button.clicked.connect(slot)

    def _get_selected_row_indices(self):
        return list(set(index.row() for index in self.grouping_table.selectedIndexes()))

    def get_selected_group_names(self):
        indexes = self._get_selected_row_indices()
        return [str(self.grouping_table.item(i, 0).getText()) for i in indexes]

    def remove_selected_groups(self):
        indices = self._get_selected_row_indices()
        for index in sorted(indices):
            self.grouping_table.removeRow(index)

    def remove_last_row(self):
        last_row = self.grouping_table.rowCount() - 1
        self.grouping_table.removeRow(last_row)

    def num_rows(self):
        return self.grouping_table.rowCount()

    def warning_popup(self, message):
        message_box.warning(str(message))
