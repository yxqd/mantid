from __future__ import (absolute_import, division, print_function)

from PyQt4 import QtCore, QtGui


class LoadView(QtGui.QWidget):
    """
    The view for the Load widget. This
    creates the look of the widget
    """
    # signals
    #helpSignal = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        super(LoadView, self).__init__(parent)
   
        self.button = QtGui.QPushButton('Load', self)
        self.button.setStyleSheet("background-color:lightgrey")

        self.button2 = QtGui.QPushButton('some data', self)
        self.button2.setStyleSheet("background-color:lightgrey")
        
        #connects
        #self.button.clicked.connect(self.helpButtonClick)
        # button layout
        self.buttonLayout=QtGui.QHBoxLayout(self)
        self.buttonLayout.addWidget(self.button)
        self.buttonLayout.addWidget(self.button2)


