from __future__ import (absolute_import, division, print_function)

import sys

import PyQt4.QtGui as QtGui
import PyQt4.QtCore as QtCore

from Muon import tab_view
from Muon import help_view
from Muon import load_view

class MainWindowView(QtGui.QWidget):
    #def __init__(self,groupedViews,groupedModels,parent=None):
    def __init__(self,parent=None):
        super(MainWindowView,self).__init__(parent)

        self.tabView = tab_view.TabView(parent)
        loadView = load_view.LoadView(parent)
        helpView = help_view.HelpView(parent)
         
        self.widget = QtGui.QWidget()
        splitter=QtGui.QSplitter(QtCore.Qt.Vertical)
        splitter.addWidget(loadView)
        splitter.addWidget(self.tabView)
        splitter.addWidget(helpView)

        splitter.setSizes([200,800,100])
        QHbox = QtGui.QHBoxLayout()
        QHbox.addWidget(splitter)
        self.widget.setLayout(QHbox)      
 
    
    def getWidget(self):
        return self.widget
   
    def closeEvent(self,event):
        self.tabView.closeEvent(event)

