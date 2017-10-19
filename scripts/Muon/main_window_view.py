from __future__ import (absolute_import, division, print_function)

import sys

import PyQt4.QtGui as QtGui
import PyQt4.QtCore as QtCore

from Muon import model_constructor
from Muon import view_constructor

from Muon import transform_presenter
from Muon import transform_view

from Muon import tab_presenter
from Muon import tab_view
from Muon import help_view
from Muon import load_view

class MainWindowView(QtGui.QWidget):
    #def __init__(self,groupedViews,groupedModels,parent=None):
    def __init__(self,parent=None):
        super(MainWindowView,self).__init__(parent)


        groupedViews = view_constructor.ViewConstructor(True,parent)
        groupedModels = model_constructor.ModelConstructor(True)


        #tabView =groupedViews.getMainView("Tabs")
        #transformView = groupedViews.getTabView("Transform")
        tabView = tab_view.TabView(parent)
        transformView = transform_view.TransformView(groupedViews,parent)
        self.presenter =tab_presenter.TabPresenter(tabView,transformView,groupedModels)
        

        #helpView = groupedViews.getMainView("Help")
        #loadView = groupedViews.getMainView("Load")
        loadView = load_view.LoadView(parent)
        helpView = help_view.HelpView(parent)
         
        self.widget = QtGui.QWidget()
        splitter=QtGui.QSplitter(QtCore.Qt.Vertical)
        splitter.addWidget(loadView)
        splitter.addWidget(tabView)
        splitter.addWidget(helpView)
        QHbox = QtGui.QHBoxLayout()
        QHbox.addWidget(splitter)
        self.widget.setLayout(QHbox)      
 
        #aself.setWindowTitle("Frequency Domain Analysis")

    def getWidget(self):
        return self.widget


