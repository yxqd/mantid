from __future__ import (absolute_import, division, print_function)


from qtpy import QtWidgets, QtGui


class BasicLoadView(QtWidgets.QWidget):

    def __init__(self, name="", parent=None):
        super(BasicLoadView, self).__init__(parent)
        self._grid = QtGui.QGridLayout(self)

        self._instrument = QtGui.QLabel(name)

        self._grid.addWidget(self.self._instrument,0,0)
        
        self._run = QtWidgets.QLineEdit()
        self._grid.addWidget(self.self._run,0,1)

        self._mode = QtWidgets.QComboBox()
        options = ["Co-add","Simultaneous"]
        self._mode.addItems(options)
        self._grid.addWidget(self.self._mode,1,0)



    def getLayout(self):
        return self.grid

