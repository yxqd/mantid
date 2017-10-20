from __future__ import (absolute_import, division, print_function)

import sys

import PyQt4.QtGui as QtGui
import PyQt4.QtCore as QtCore

import collections

from Muon import tab_view
from Muon import help_view
from Muon import load_view

class MainWindowView(QtGui.QWidget):
    def __init__(self,parent=None):
        super(MainWindowView,self).__init__(parent)

        self.views = self.populateMainWidget(parent)
         
        self.widget = QtGui.QWidget()
        splitter=QtGui.QSplitter(QtCore.Qt.Vertical)
        
        for key in self.views:
             splitter.addWidget(self.views[key])

        splitter.setSizes(self.getSizes())
        QHbox = QtGui.QHBoxLayout()
        QHbox.addWidget(splitter)
        self.widget.setLayout(QHbox)      
 
    
    def getWidget(self):
        return self.widget
 
    def getViews(self):
        return self.views
  
    def closeEvent(self,event):
        for key in self.views:
            self.views[key].closeEvent(event)
   
    """
    Only the following code would need updating for a new GUI
    """
    def populateMainWidget(self,parent):
        """
        This code may need updating for a different GUI
        It constructs the different views that are contained 
        within the main window
   
        The order must match the order they appear in the window
        (top to bottom)
        """
        views = collections.OrderedDict()
        views["Load"]=(load_view.LoadView(parent))
        views["Tab"]=(tab_view.TabView(parent))
        views["Help"]=(help_view.HelpView(parent))
        return views 

    def getSizes(self):
        return [200,800,100]
