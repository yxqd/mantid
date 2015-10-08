#pylint: disable=invalid-name

from PyQt4 import QtCore, QtGui
from ui_ansto_sans import Ui_WindowSansAnsto
from view.main_view import MainView
from view.processing_view import ProcessingView
from view.settings_view import SettingsView
from view.masking_view import MaskingView

from presenter.main_view_presenter import MainViewPresenter

class QtProcessingView(ProcessingView):

    def __init__(self, processingTab):
        self._processingTab = processingTab

    def show_empty(self):
        '''
        layout = QtGui.QVBoxLayout()
        layout.addWidget(QtGui.QLineEdit())
        self._processingTab.setLayout(layout)
        '''
        pass

class QtSettingsView(SettingsView):

    def __init__(self, settingsTab):
        self._settingsTab = settingsTab

    def show_empty(self):
        pass

class QtMaskingView(MaskingView):
    def __init__(self, settingsTab):
        self._settingsTab = settingsTab

    def show_empty(self):
        pass





class SANSBatchGui(MainView, QtGui.QMainWindow, Ui_WindowSansAnsto):


    def __init__(self):
        """
            Constructor
        """
        super(QtGui.QMainWindow, self).__init__()
        self.setupUi(self)

        _presenter = MainViewPresenter(self)


    def setup_layout(self):
        return True

    def select_processing_view(self):
        self.tabs.setCurrentIndex(0)

    def get_processing_view(self):
        return QtProcessingView(self.processing)

    def get_settings_view(self):
        return QtSettingsView(self.settings)

    def get_mask_view(self):
        return QtMaskingView(self.masking)
