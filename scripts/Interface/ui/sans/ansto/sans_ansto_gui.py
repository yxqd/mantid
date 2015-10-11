#pylint: disable=invalid-name

from PyQt4 import QtGui

from ui_ansto_sans import Ui_WindowSansAnsto
from view.main_view import MainView
from scripts.Interface.ui.sans.ansto.main_view_presenter import MainViewPresenter
from command import Command


class SANSBatchGui(MainView, QtGui.QMainWindow, Ui_WindowSansAnsto):


    def __init__(self):
        """
            Constructor
        """
        super(QtGui.QMainWindow, self).__init__()
        self.setupUi(self)

        self._presenter = MainViewPresenter(self)


        self.btn_reduce.clicked.connect(self.on_reduce)

    def on_reduce(self):
        self._presenter.notify(Command.ProcessAll)

    def setup_layout(self):
        return True

    def set_processing(self):
        print "Processing!"


