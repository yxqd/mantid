#pylint: disable=invalid-name

from PyQt4 import QtGui

from ui_ansto_sans import Ui_WindowSansAnsto
from view.main_view import MainView
from main_view_presenter import MainViewPresenter
from command import Command

class SansBatchGuiMetaClass(type(MainView), type(QtGui.QMainWindow)):
    pass

class SANSBatchGui(MainView, QtGui.QMainWindow, Ui_WindowSansAnsto):

    __metaclass__ = SansBatchGuiMetaClass

    def __init__(self):
        """
            Constructor
        """
        super(QtGui.QMainWindow, self).__init__()
        self.setupUi(self)
        # Create the presenter/controller for the UI
        self._presenter = MainViewPresenter(self)
        # For the click processing event
        self.btn_reduce.clicked.connect(self.on_reduce)

    def on_reduce(self):
        # Update the presenter.
        self._presenter.notify(Command.ProcessAll)

    def setup_layout(self):
        return True

    def set_processing(self):
        # TODO modal dialogs etc
        print "Processing Runs"

    def get_run_table(self):
        # TODO Get all run information
        return None


