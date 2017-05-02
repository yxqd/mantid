from __future__ import (absolute_import, division, print_function)
import Indirect_QuickRun.q
from QuickRun.quickrun_presenter import QuickRunPresenter
from QuickRun.quickrun_model import PlotOptionsModel
from PyQt4.QtGui import QApplication
import sys


def main():
    app = QApplication(sys.argv)
    view = QuickRunView()
    plotModel = PlotOptionsModel()
    presenter = QuickRunPresenter(view, plotModel)
    view.show()
    app.exec_()

