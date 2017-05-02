from PyQt4 import uic, QtCore
from PyQt4.QtGui import QMainWindow, QApplication, QWidget
import sys
from quickrun_presenter import QuickRunPresenter
from quickrun_model import PlotOptionsModel


class ScanTab(QWidget):
    def __init__(self):
        super(ScanTab, self).__init__()
        self.energytab = uic.loadUi('./EnergyWindowScan.ui')
        self.sqwtab = uic.loadUi('./SQWMomentScan.ui')
        self.diffractiontab = uic.loadUi('./DiffractionScan.ui')
        self.samplechangertab = uic.loadUi('./SampleChanger.ui')


class QuickRunView(QMainWindow):
    plotClicked = QtCore.pyqtSignal()
    tabChanged = QtCore.pyqtSignal(int)

    def __init__(self):
        super(MainWindowView, self).__init__()
        self.ui = uic.loadUi('./IndirectQuickRun.ui', self)
        self.setWindowTitle("Indirect QuickRun")
        self.tabs = ScanTab()
        self.ui.tb_quickrun.insertTab(0, self.tabs.energytab, "Energy Window Scan")
        self.ui.tb_quickrun.insertTab(1, self.tabs.sqwtab, "SQW Moments Scan")
        self.ui.tb_quickrun.insertTab(2, self.tabs.diffractiontab, "Diffraction Scan")
        self.ui.tb_quickrun.insertTab(3, self.tabs.samplechangertab, "Sample Changer")

        # set default plot options
        self.ui.cb_plotOptions.addItems(['Spectra', 'Contour', 'Elwin', 'MSDFit'])

        # connect slots
        self.ui.pb_plot.clicked.connect(self._onPlot)
        self.ui.tb_quickrun.currentChanged.connect(self._onTab)

    def _onPlot(self):
        self.plotClicked.emit()

    def _onTab(self, index):
        self.tabChanged.emit(index)

    def getPlotOption(self):
        """
        returns the currently selected plot option
        """
        return self.ui.cb_plotOptions.currentText()

    def addPlotOptions(self, items):
        """
        sets the plot option combo box to a list of options
        :param items: list of plot options
        """

        self.ui.cb_plotOptions.clear()
        self.ui.cb_plotOptions.addItems(items)
        if not items:
            self.ui.cb_plotOptions.setDisabled(True)
        else:
            self.ui.cb_plotOptions.setDisabled(False)
        return None
