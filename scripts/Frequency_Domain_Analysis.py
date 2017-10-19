#pylint: disable=invalid-name
from __future__ import (absolute_import, division, print_function)

import sys

import PyQt4.QtGui as QtGui
import PyQt4.QtCore as QtCore

from Muon import model_constructor
from Muon import view_constructor
from Muon import main_window_view

class FrequencyDomainAnalysisGui(QtGui.QMainWindow):
    def __init__(self,parent=None):
        super(FrequencyDomainAnalysisGui,self).__init__(parent)

        #groupedViews = view_constructor.ViewConstructor(True,self)
        #groupedModels = model_constructor.ModelConstructor(True)
        widget = main_window_view.MainWindowView(parent)
        self.setCentralWidget(widget.getWidget())
        self.setWindowTitle("Frequency Domain Analysis")

    # cancel algs if window is closed
    def closeEvent(self,event):
        a=1
        #self.presenter.close()

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
