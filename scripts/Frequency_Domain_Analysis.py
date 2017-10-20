#pylint: disable=invalid-name
from __future__ import (absolute_import, division, print_function)

import sys

import PyQt4.QtGui as QtGui
import PyQt4.QtCore as QtCore

from Muon import main_window_view
from Muon import main_window_presenter

from Muon import model_constructor

class FrequencyDomainAnalysisGui(QtGui.QMainWindow):
    def __init__(self,parent=None):
        super(FrequencyDomainAnalysisGui,self).__init__(parent)

        groupedModels = model_constructor.ModelConstructor(True)

        self.view = main_window_view.MainWindowView(self)
        self.presenter=main_window_presenter.MainWindowPresenter(self.view,groupedModels)
        self.setCentralWidget(self.view.getWidget())
        self.setWindowTitle("Frequency Domain Analysis")

    # cancel algs if window is closed
    def closeEvent(self,event):
        self.presenter.close()
        self.view.closeEvent(event)

def qapp():
    if QtGui.QApplication.instance():
        _app = QtGui.QApplication.instance()
    else:
        _app = QtGui.QApplication(sys.argv)
    return _app


app = qapp()
ex= FrequencyDomainAnalysisGui()
ex.resize(700,700)
ex.show()
app.exec_()
