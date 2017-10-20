from __future__ import (absolute_import, division, print_function)

from PyQt4 import QtGui
from PyQt4 import QtCore


from Muon import view_constructor
from Muon import transform_view

class TabView(QtGui.QMainWindow):
    """
    Creates the view for the tabs.
    """
    def __init__(self,parent=None):
        super(TabView,self).__init__(parent)
        #self.dock = QtGui.QMainWindow()
        self.widgets=[]
        self.tabs=[]
        groupedViews=view_constructor.ViewConstructor(True,self)
        self.transformView = transform_view.TransformView(groupedViews)

    def addTab(self,widget,name):
        # add widget
        self.widgets.append(widget)
        # make an empty dock for widget
        self.tabs.append(QtGui.QDockWidget(name))
        # add widget to dock
        self.tabs[-1].setWidget(self.widgets[-1])
        #add dock to view
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea,self.tabs[-1])
        # set layout
        self.widgets[-1].setLayout(self.widgets[-1].getLayout())

    def makeTabs(self):
        # convert docks into tabs
        self.tabifyDockWidget(self.tabs[0],self.tabs[1])
        for j in range(2,len(self.tabs),1):
            self.tabifyDockWidget(self.tabs[0],self.tabs[j])
        # put tabs on the top of the page 
        self.setTabPosition(QtCore.Qt.LeftDockWidgetArea,QtGui.QTabWidget.North)
        # open to the first tab
        self.tabs[0].show()
        self.tabs[0].raise_()

        # stop tabs from closing
        for j in range(0,len(self.tabs),1):
           self.tabs[j].setFeatures(QtGui.QDockWidget.DockWidgetClosable and QtGui.QDockWidget.DockWidgetFloatable)

    def closeEvent(self,event):
        for j in range(len(self.tabs)-1,-1,-1):
            self.tabs[j].close()
            self.widgets[j].close()
