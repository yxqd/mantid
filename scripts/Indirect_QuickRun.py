from __future__ import (absolute_import, division, print_function)
from IndirectQuickRun.quickrun_view import QuickRunView
from IndirectQuickRun.quickrun_presenter import QuickRunPresenter
from IndirectQuickRun.quickrun_model import PlotOptionsModel
from PyQt4.QtGui import QApplication
import sys


def main():
    app = QApplication(sys.argv)
    view = QuickRunView()
    plotModel = PlotOptionsModel()
    presenter = QuickRunPresenter(view, plotModel)
    view.show()
    app.exec_()

if __name__ == '__main__':
    main()