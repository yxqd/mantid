from __future__ import (absolute_import, division, print_function)
import sys
from PyQt4 import QtGui

from Muon.GUI.Common.grouping_table_widget.grouping_table_widget_model import GroupingTableModel
from Muon.GUI.Common.grouping_table_widget.grouping_table_widget_view import GroupingTableView
from Muon.GUI.Common.grouping_table_widget.grouping_table_widget_presenter import GroupingTablePresenter, MuonGroup

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    ui = GroupingTablePresenter(GroupingTableView(), GroupingTableModel())

    testgroup1 = MuonGroup(group_name="fwd", detector_IDs=[1, 2, 3, 4, 5])
    testgroup2 = MuonGroup(group_name="bwd", detector_IDs=[6, 7, 8, 9, 10])
    testgroup3 = MuonGroup(group_name="top", detector_IDs=[11, 12, 13, 14, 15])
    ui.add_group(testgroup1)
    ui.add_group(testgroup2)
    ui.add_group(testgroup3)

    ui.show()
    sys.exit(app.exec_())
