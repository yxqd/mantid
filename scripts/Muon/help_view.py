from __future__ import (absolute_import, division, print_function)

from PyQt4 import QtCore, QtGui


class HelpView(QtGui.QWidget):
    """
    The view for the Help widget. This
    creates the look of the widget
    """
    # signals
    #helpSignal = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        super(HelpView, self).__init__(parent)
  
        self.button = QtGui.QPushButton('?', self)
        self.button.setStyleSheet("background-color:lightgrey")
        #self.spacer= QtGui.QSpacerItem(20, 30, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        #connects
        #self.button.clicked.connect(self.helpButtonClick)
        # button layout
        self.buttonLayout=QtGui.QGridLayout(self)#HBoxLayout(self)
        self.buttonLayout.addWidget(self.button)
        #self.buttonLayout.addWidget(self.spacer)

    def getLayout(self):
        return self.buttonLayout
