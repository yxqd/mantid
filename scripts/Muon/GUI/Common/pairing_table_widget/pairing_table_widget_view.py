from __future__ import (absolute_import, division, print_function)

from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import pyqtSignal as Signal

from functools import wraps

from Muon.GUI.Common import message_box


class PairNameTableItem(QtGui.QTableWidgetItem):
    """
    An extension of the QTableWidgetItem class, which modifies the setData method to first check that the enetred
    text is valid; and only runs setData if the validator returns True
    """

    @staticmethod
    def validator_before_set(func, validator):
        @wraps(func)
        def wrapper(*args, **kw):
            if validator(args[1].toString()):
                res = func(*args, **kw)
            else:
                res = None
            return res

        return wrapper

    @staticmethod
    def default_validator(text):
        return True

    def __init__(self, validator=default_validator):
        """
        :param validator: A predicate fucntion (returns True/False) taking a single string as argument
        """
        super(PairNameTableItem, self).__init__(0)
        self._validator = validator
        self._modify_setData()

    def validator(self, text):
        return self._validator(text)

    def _modify_setData(self):
        """Modify the setData method"""
        setattr(self, "setData", self.validator_before_set(self.setData, self.validator))


class PairingTableView(QtGui.QWidget):
    dataChanged = Signal()

    def __init__(self, parent=None):
        super(PairingTableView, self).__init__(parent)

        self.pairing_table = QtGui.QTableWidget(self)
        self.set_up_table()

        self.setup_interface_layout()

        self.pairing_table.itemChanged.connect(self.on_item_changed)
        self.pairing_table.cellChanged.connect(self.on_cell_changed)

        self._validate_pair_name_entry = lambda text: True

        self._on_table_data_changed = lambda: 0

        self.cached_table = None

        # The active groups that can be selected from the group combo box
        self._group_selections = []

        # whether the table is updating and therefore
        # we shouldn't respond to signals
        self._updating = False

    def update_group_selections(self, group_name_list):
        self._group_selections = group_name_list

    def on_user_changes_pair_name(self, slot):
        self._validate_pair_name_entry = slot

    def setup_interface_layout(self):
        self.setObjectName("PairingTableView")
        self.resize(500, 500)

        self.add_pair_button = QtGui.QToolButton()
        self.remove_pair_button = QtGui.QToolButton()

        size_policy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        size_policy.setHorizontalStretch(0)
        size_policy.setVerticalStretch(0)
        size_policy.setHeightForWidth(self.add_pair_button.sizePolicy().hasHeightForWidth())
        size_policy.setHeightForWidth(self.remove_pair_button.sizePolicy().hasHeightForWidth())

        self.add_pair_button.setSizePolicy(size_policy)
        self.add_pair_button.setMinimumSize(QtCore.QSize(25, 25))
        self.add_pair_button.setObjectName("addGroupButton")
        self.add_pair_button.setText("+")

        self.remove_pair_button.setSizePolicy(size_policy)
        self.remove_pair_button.setMinimumSize(QtCore.QSize(25, 25))
        self.remove_pair_button.setObjectName("removeGroupButton")
        self.remove_pair_button.setText("-")

        self.horizontal_layout = QtGui.QHBoxLayout()
        self.horizontal_layout.setObjectName("horizontalLayout")
        self.horizontal_layout.addWidget(self.add_pair_button)
        self.horizontal_layout.addWidget(self.remove_pair_button)
        self.spacer_item = QtGui.QSpacerItem(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.horizontal_layout.addItem(self.spacer_item)
        self.horizontal_layout.setAlignment(QtCore.Qt.AlignLeft)

        self.vertical_layout = QtGui.QVBoxLayout(self)
        self.vertical_layout.setObjectName("verticalLayout")
        self.vertical_layout.addWidget(self.pairing_table)
        self.vertical_layout.addLayout(self.horizontal_layout)

        self.setLayout(self.vertical_layout)

    def cache_table(self):
        self.cached_table = self.get_table_contents()

    def restore_cached_state(self):
        cache = self.cached_table
        if cache:
            for i in range(len(cache)):
                pass

    def set_up_table(self):
        self.pairing_table.setColumnCount(4)
        self.pairing_table.setHorizontalHeaderLabels(QtCore.QString("Pair Name;Group 1; Group 2;Alpha").split(";"))
        header = self.pairing_table.horizontalHeader()
        header.setResizeMode(0, QtGui.QHeaderView.Stretch)
        header.setResizeMode(1, QtGui.QHeaderView.Stretch)
        header.setResizeMode(2, QtGui.QHeaderView.Stretch)
        header.setResizeMode(3, QtGui.QHeaderView.ResizeToContents)
        # table_utils.setTableHeaders(self.grouping_table)
        vertical_headers = self.pairing_table.verticalHeader()
        vertical_headers.setMovable(False)
        vertical_headers.setResizeMode(QtGui.QHeaderView.ResizeToContents)
        vertical_headers.setVisible(True)

    def contextMenuEvent(self, event):
        self.menu = QtGui.QMenu(self)
        add_pair_action = QtGui.QAction('Add Pair', self)
        add_pair_action.triggered.connect(self.add_pair_button.clicked.emit)

        remove_pair_action = QtGui.QAction('Remove Pair', self)
        remove_pair_action.triggered.connect(self.remove_pair_button.clicked.emit)

        self.menu.addAction(add_pair_action)
        self.menu.addAction(remove_pair_action)
        # add other required actions
        self.menu.popup(QtGui.QCursor.pos())

    def _group_selection_cell_widget(self):
        selector = QtGui.QComboBox(self)
        selector.addItems(self._group_selections)
        # selector.setStyleSheet("border: 0px solid white;")
        return selector

    def get_index_of_text(self, selector, text):
        for i in range(selector.count()):
            if str(selector.itemText(i)) == text:
                return i
        return 0

    def add_entry_to_table(self, row_entries):
        assert len(row_entries) == self.pairing_table.columnCount()

        row_position = self.pairing_table.rowCount()
        self.pairing_table.insertRow(row_position)
        for i, entry in enumerate(row_entries):
            item = QtGui.QTableWidgetItem(entry)
            if i == 0:
                pair_name_widget = PairNameTableItem(self._validate_pair_name_entry)
                pair_name_widget.setText(entry)
                self.pairing_table.setItem(row_position, 0, pair_name_widget)
                continue
            if i == 1:
                group1_selector_widget = self._group_selection_cell_widget()
                index = self.get_index_of_text(group1_selector_widget, entry)
                group1_selector_widget.setCurrentIndex(index)
                self.pairing_table.setCellWidget(row_position, 1, group1_selector_widget)
                continue
            if i == 2:
                group2_selector_widget = self._group_selection_cell_widget()
                index = self.get_index_of_text(group2_selector_widget, entry)
                group2_selector_widget.setCurrentIndex(index)
                self.pairing_table.setCellWidget(row_position, 2, group2_selector_widget)
                continue
            if i == 3:
                item.setFlags(QtCore.Qt.ItemIsEnabled)
                item.setFlags(QtCore.Qt.ItemIsSelectable)
            self.pairing_table.setItem(row_position, i, item)

    def on_add_pair_button_clicked(self, slot):
        self.add_pair_button.clicked.connect(slot)

    def on_remove_pair_button_clicked(self, slot):
        self.remove_pair_button.clicked.connect(slot)

    def on_table_data_changed(self, slot):
        self._on_table_data_changed = slot

    def _get_selected_row_indices(self):
        return list(set(index.row() for index in self.pairing_table.selectedIndexes()))

    def get_selected_pair_names(self):
        indexes = self._get_selected_row_indices()
        return [str(self.pairing_table.item(i, 0).text()) for i in indexes]

    def remove_selected_pairs(self):
        indices = self._get_selected_row_indices()
        for index in reversed(sorted(indices)):
            self.pairing_table.removeRow(index)

    def remove_last_row(self):
        last_row = self.pairing_table.rowCount() - 1
        self.pairing_table.removeRow(last_row)

    def num_rows(self):
        return self.pairing_table.rowCount()

    def warning_popup(self, message):
        message_box.warning(str(message))

    def on_item_changed(self):
        if not self._updating:
            return
        return

    def on_cell_changed(self, row, col):
        if not self._updating:
            self._on_table_data_changed()

    def get_table_contents(self):
        if self._updating:
            return []
        ret = [[None for _ in range(4)] for _ in range(self.num_rows())]
        for i in range(self.num_rows()):
            for j in range(4):
                if j == 1 or j == 2:
                    ret[i][j] = str(self.pairing_table.cellWidget(i, j).currentText())
                else:
                    ret[i][j] = str(self.pairing_table.item(i, j).text())

        return ret

    def clear(self):
        for row in reversed(range(self.num_rows())):
            self.pairing_table.removeRow(row)

    def notify_data_changed(self):
        self.dataChanged.emit()

    def disable_updates(self):
        self._updating = True

    def enable_updates(self):
        self._updating = False
