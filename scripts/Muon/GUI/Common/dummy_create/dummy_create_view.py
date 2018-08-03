from __future__ import (absolute_import, division, print_function)

from PyQt4 import QtCore, QtGui

import mantid.simpleapi as mantid

from Muon.GUI.Common import table_utils


class DummyCreateView(QtGui.QWidget):

    """
    creates the layout for the FFT GUI
    """
    # signals
    buttonSignal = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        super(DummyCreateView, self).__init__(parent)
        self.grid = QtGui.QGridLayout(self)

        # add splitter for resizing
        splitter = QtGui.QSplitter(QtCore.Qt.Vertical)

        # make table
        self.Table = QtGui.QTableWidget(self)
        self.Table.resize(800, 800)
        self.Table.setRowCount(3)
        self.Table.setColumnCount(2)
        self.Table.setColumnWidth(0, 300)
        self.Table.setColumnWidth(1, 300)
        self.Table.verticalHeader().setVisible(False)
        self.Table.horizontalHeader().setStretchLastSection(True)
        self.Table.setHorizontalHeaderLabels(
            ("Property;Value").split(";"))
        # populate table
        self.widgets = {}

        table_utils.setRowName(self.Table, 0, "x data")
        self.widgets["DataX"] = table_utils.addStringToTable(self.Table, "0.0,1.0,2.0", 0)
  
        table_utils.setRowName(self.Table, 1, "y data")
        self.widgets["DataY"] = table_utils.addStringToTable(self.Table, "0.0,-1.0,-2.0", 1)

        table_utils.setRowName(self.Table, 2, "workspace name")
        self.widgets["OutputWorkspace"] = table_utils.addStringToTable(self.Table, "Muon", 2)

        # make button
        self.button = QtGui.QPushButton('Create', self)
        self.button.setStyleSheet("background-color:lightgrey")
        # connects
        self.button.clicked.connect(self.buttonClick)
        # add to layout
        self.grid.addWidget(self.Table)
        self.grid.addWidget(self.button)

    def getLayout(self):
        return self.grid

    def buttonClick(self):
        self.buttonSignal.emit()

    def getDataX(self):
        return self.widgets["DataX"].text()

    def getDataY(self):
        return self.widgets["DataY"].text()

    def getWS(self):
        return self.widgets["OutputWorkspace"].text()
