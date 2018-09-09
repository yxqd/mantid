from __future__ import (absolute_import, division, print_function)

from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import pyqtSignal as Signal

from functools import wraps
import re

from Muon.GUI.Common import table_utils
from Muon.GUI.Common import message_box


def decorate_setData_with_regex(function_decorator, regex):
    def decorator(cls):
        if "setData" in dir(cls):
            obj = getattr(cls, "setData")
            print("Decorator obj is : ", obj)
            try:
                obj = obj.__func__  # unwrap Python 2 unbound method
            except AttributeError:
                pass  # not needed in Python 3
            setattr(cls, "setData", function_decorator(obj, regex))
        return cls

    return decorator


def regex_before_set(func, regex):
    @wraps(func)
    def wrapper(*args, **kw):
        # print('{} called'.format(func.__name__))
        # print(args[2].toString())
        # print(bool(re.match(regex, args[2].toString())))
        if bool(re.match(regex, args[2].toString())):
            res = func(*args, **kw)
        else:
            res = None
        return res

    return wrapper


class GroupNameTableItem(QtGui.QTableWidgetItem):
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
        super(GroupNameTableItem, self).__init__(0)
        self._validator = validator
        self._modify_setData()

    def validator(self, text):
        return self._validator(text)

    def _modify_setData(self):
        """Modify the setData method"""
        setattr(self, "setData", self.validator_before_set(self.setData, self.validator))


@decorate_setData_with_regex(regex_before_set, "^[0-9]*([0-9]+[,-]{0,1})*[0-9]+$")
class DetectorTableItem(QtGui.QTableWidgetItem):

    def __init__(self, validator=lambda text: True):
        super(DetectorTableItem, self).__init__(0)

        self.validator = validator


class GroupingTableView(QtGui.QWidget):
    dataChanged = Signal()

    def __init__(self, parent=None):
        super(GroupingTableView, self).__init__(parent)

        self.grouping_table = QtGui.QTableWidget(self)
        self.set_up_table()

        self.setup_interface_layout()

        self.grouping_table.itemChanged.connect(self.on_item_changed)
        self.grouping_table.cellChanged.connect(self.on_cell_changed)

        self._validate_group_name_entry = lambda text: True
        self._validate_detector_ID_entry = lambda text: True

        self._on_table_data_changed = lambda: 0

        self.cached_table = None

        # whether the table is updating and therefore
        # we shouldn't respond to signals
        self._updating = False

    def on_user_changes_group_name(self, slot):
        self._validate_group_name_entry = slot

    def on_user_changes_detector_IDs(self, slot):
        self._validate_detector_ID_entry = slot

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

        self.horizontal_layout = QtGui.QHBoxLayout()
        self.horizontal_layout.setObjectName("horizontalLayout")
        self.horizontal_layout.addWidget(self.add_group_button)
        self.horizontal_layout.addWidget(self.remove_group_button)
        self.spacer_item = QtGui.QSpacerItem(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.horizontal_layout.addItem(self.spacer_item)
        self.horizontal_layout.setAlignment(QtCore.Qt.AlignLeft)

        self.vertical_layout = QtGui.QVBoxLayout(self)
        self.vertical_layout.setObjectName("verticalLayout")
        self.vertical_layout.addWidget(self.grouping_table)
        self.vertical_layout.addLayout(self.horizontal_layout)
        # self.vertical_layout.addWidget(self.add_group_button)
        # self.vertical_layout.addWidget(self.remove_group_button)

        self.setLayout(self.vertical_layout)

    def cache_table(self):
        self.cached_table = self.get_table_contents()

    def restore_cached_state(self):
        cache = self.cached_table
        if cache:
            for i in range(len(cache)):
                pass

    def set_up_table(self):
        self.grouping_table.setColumnCount(3)
        self.grouping_table.setHorizontalHeaderLabels(QtCore.QString("Group Name;Detector IDs;N Detectors").split(";"))
        header = self.grouping_table.horizontalHeader()
        header.setResizeMode(0, QtGui.QHeaderView.Stretch)
        header.setResizeMode(1, QtGui.QHeaderView.Stretch)
        header.setResizeMode(2, QtGui.QHeaderView.ResizeToContents)
        # table_utils.setTableHeaders(self.grouping_table)
        vertical_headers = self.grouping_table.verticalHeader()
        vertical_headers.setMovable(False)
        vertical_headers.setResizeMode(QtGui.QHeaderView.ResizeToContents)
        vertical_headers.setVisible(True)

    def contextMenuEvent(self, event):
        self.menu = QtGui.QMenu(self)
        addGroupAction = QtGui.QAction('Add Group', self)
        addGroupAction.triggered.connect(self.add_group_button.clicked.emit)

        removeGroupAction = QtGui.QAction('Remove Group', self)
        removeGroupAction.triggered.connect(self.remove_group_button.clicked.emit)

        self.menu.addAction(addGroupAction)
        self.menu.addAction(removeGroupAction)
        # add other required actions
        self.menu.popup(QtGui.QCursor.pos())

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
                group_name_widget = GroupNameTableItem(self._validate_group_name_entry)
                group_name_widget.setText(entry)
                self.grouping_table.setItem(row_position, 0, group_name_widget)
                continue
            if i == 1:
                detector_widget = DetectorTableItem()  # self._detector_cell_widget()
                detector_widget.setText(entry)
                self.grouping_table.setItem(row_position, 1, detector_widget)
                continue
            if i == 2:
                item.setFlags(QtCore.Qt.ItemIsEnabled)
                item.setFlags(QtCore.Qt.ItemIsSelectable)
            self.grouping_table.setItem(row_position, i, item)

    def on_add_group_button_clicked(self, slot):
        self.add_group_button.clicked.connect(slot)

    def on_remove_group_button_clicked(self, slot):
        self.remove_group_button.clicked.connect(slot)

    def on_table_data_changed(self, slot):
        self._on_table_data_changed = slot

    def _get_selected_row_indices(self):
        return list(set(index.row() for index in self.grouping_table.selectedIndexes()))

    def get_selected_group_names(self):
        indexes = self._get_selected_row_indices()
        return [str(self.grouping_table.item(i, 0).text()) for i in indexes]

    def remove_selected_groups(self):
        indices = self._get_selected_row_indices()
        for index in reversed(sorted(indices)):
            self.grouping_table.removeRow(index)

    def remove_last_row(self):
        last_row = self.grouping_table.rowCount() - 1
        self.grouping_table.removeRow(last_row)

    def num_rows(self):
        return self.grouping_table.rowCount()

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
        ret = [[None for _ in range(3)] for _ in range(self.num_rows())]
        for i in range(self.num_rows()):
            for j in range(3):
                ret[i][j] = str(self.grouping_table.item(i, j).text())
        return ret

    def notify_data_changed(self):
        self.dataChanged.emit()

    def clear(self):
        for row in reversed(range(self.num_rows())):
            self.grouping_table.removeRow(row)

    def disable_updates(self):
        self._updating = True

    def enable_updates(self):
        self._updating = False
