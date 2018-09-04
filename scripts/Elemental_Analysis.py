from __future__ import absolute_import, print_function

from PyQt4 import QtGui

import sys
import random
from time import time

from Muon.GUI.ElementalAnalysis.PeriodicTable.periodic_table_presenter import PeriodicTablePresenter
from Muon.GUI.ElementalAnalysis.PeriodicTable.periodic_table_view import PeriodicTableView
from Muon.GUI.ElementalAnalysis.PeriodicTable.periodic_table_model import PeriodicTableModel
from Muon.GUI.ElementalAnalysis.Plotting.plotting_view import PlotView
from Muon.GUI.ElementalAnalysis.Plotting.plotting_presenter import PlotPresenter
from Muon.GUI.Common import message_box

from Muon.GUI.ElementalAnalysis.Detectors.detectors_presenter import DetectorsPresenter
from Muon.GUI.ElementalAnalysis.Detectors.detectors_view import DetectorsView
from Muon.GUI.ElementalAnalysis.Peaks.peaks_presenter import PeaksPresenter
from Muon.GUI.ElementalAnalysis.Peaks.peaks_view import PeaksView

from Muon.GUI.ElementalAnalysis.PeriodicTable.PeakSelector.peak_selector_presenter import PeakSelectorPresenter
from Muon.GUI.ElementalAnalysis.PeriodicTable.PeakSelector.peak_selector_view import PeakSelectorView
from Muon.GUI.ElementalAnalysis.LoadWidget.load_model import LoadModel, CoLoadModel
from Muon.GUI.Common.load_widget.load_view import LoadView
from Muon.GUI.Common.load_widget.load_presenter import LoadPresenter

import mantid.simpleapi as mantid

import matplotlib.patches as mpatches


class ElementalAnalysisGui(QtGui.QMainWindow):
    def __init__(self, parent=None):
        super(ElementalAnalysisGui, self).__init__(parent)
        mantid.LoadAscii(
            "~/Sundials/ral02695.rooth2010.dat",
            OutputWorkspace="D1N10")
        mantid.LoadAscii(
            "~/Sundials/ral02695.rooth3010.dat",
            OutputWorkspace="D2N10")
        mantid.LoadAscii(
            "~/Sundials/ral02695.rooth4010.dat",
            OutputWorkspace="D3N10")
        mantid.LoadAscii(
            "~/Sundials/ral02695.rooth5010.dat",
            OutputWorkspace="D4N10")

        mantid.LoadAscii(
            "~/Sundials/ral02695.rooth2020.dat",
            OutputWorkspace="D1N20")
        mantid.LoadAscii(
            "~/Sundials/ral02695.rooth3020.dat",
            OutputWorkspace="D2N20")
        mantid.LoadAscii(
            "~/Sundials/ral02695.rooth4020.dat",
            OutputWorkspace="D3N20")
        mantid.LoadAscii(
            "~/Sundials/ral02695.rooth5020.dat",
            OutputWorkspace="D4N20")

        mantid.LoadAscii(
            "~/Sundials/ral02695.rooth2099.dat",
            OutputWorkspace="D1N99")
        mantid.LoadAscii(
            "~/Sundials/ral02695.rooth3099.dat",
            OutputWorkspace="D2N99")
        mantid.LoadAscii(
            "~/Sundials/ral02695.rooth4099.dat",
            OutputWorkspace="D3N99")
        mantid.LoadAscii(
            "~/Sundials/ral02695.rooth5099.dat",
            OutputWorkspace="D4N99")

        self.menu = self.menuBar()
        self.menu.addAction("File")
        edit_menu = self.menu.addMenu("Edit")
        edit_menu.addAction("Change Peak Data file", self.select_data_file)
        self.menu.addAction("Binning")
        self.menu.addAction("Normalise")

        self.ptable = PeriodicTablePresenter(
            PeriodicTableView(), PeriodicTableModel())
        self.ptable.register_table_changed(self.table_changed)
        self.ptable.register_table_lclicked(self.table_left_clicked)
        self.ptable.register_table_rclicked(self.table_right_clicked)

        self.widget_list = QtGui.QVBoxLayout()

        self.load_widget = LoadPresenter(
            LoadView(), LoadModel(), CoLoadModel())
        self.plotting = PlotPresenter(PlotView())
        self.plotting.view.setMinimumSize(self.plotting.view.sizeHint())

        self.detectors = DetectorsPresenter(DetectorsView())
        for checkbox in [self.detectors.view.GE1, self.detectors.view.GE2,
                         self.detectors.view.GE3, self.detectors.view.GE4]:
            checkbox.on_checkbox_checked(self.add_checkbox_plot)
            checkbox.on_checkbox_unchecked(self.del_checkbox_plot)
        self.peaks = PeaksPresenter(PeaksView())

        self.widget_list.addWidget(self.peaks.view)
        self.widget_list.addWidget(self.detectors.view)
        self.widget_list.addWidget(self.load_widget.view)

        self.box = QtGui.QHBoxLayout()
        self.box.addWidget(self.ptable.view)
        self.box.addLayout(self.widget_list)
        self.setCentralWidget(QtGui.QWidget(self))
        self.centralWidget().setLayout(self.box)
        self.setWindowTitle("Elemental Analysis")

        self.element_widgets = {}
        self._generate_element_widgets()

        # for detector in [1, 2, 3, 4]:
        #     names = ["D{}N{}".format(detector, x) for x in [10, 20, 99]]
        #     self.plot("Detector {}".format(detector), names)
        labels = ["Zn"]  # , "Cu"]
        colours = ["b"]  # , "r"]
        mps = [mpatches.Patch(color=c, label=l)
               for c, l in zip(colours, labels)]
        self.plotting.view.figure.legend(
            mps, labels, loc="top right", prop={"size": 10})
        # self.plotting.view.figure.tight_layout()

    def add_checkbox_plot(self, checkbox):
        detector = checkbox.name[-1]
        names = ["D{}N{}".format(detector, x) for x in [10, 20, 99]]
        self.plot(checkbox.name, names)
        if self.plotting.view.isHidden():
            self.plotting.view.show()

    def del_checkbox_plot(self, checkbox):
        self.plotting.remove_subplot(checkbox.name)
        if not len(self.plotting.get_subplots()):
            self.plotting.view.close()

    def _generate_element_widgets(self):
        self.element_widgets = {}
        for element in self.ptable.peak_data:
            if element not in ["Gammas", "Electrons"]:
                data = self.ptable.element_data(element)
                widget = PeakSelectorPresenter(PeakSelectorView(data, element))
                self.element_widgets[element] = widget

    def table_left_clicked(self, item):
        print("Element Left Clicked: {}".format(
            self.ptable.element_data(item.symbol)))

    def table_right_clicked(self, item):
        self.element_widgets[item.symbol].view.show()
        print("Element Right Clicked: {}".format(item.symbol))

    def table_changed(self, items):
        print("Table Changed: {}".format([i.symbol for i in items]))

    def select_data_file(self):
        filename = str(QtGui.QFileDialog.getOpenFileName())
        if filename:
            self.ptable.set_peak_datafile(filename)

    def plot(self, detector_name, workspaces):
        subplot = self.plotting.add_subplot(detector_name)
        self.plotting.call_plot_method(
            detector_name, subplot.set_title, detector_name)
        for workspace in workspaces:
            self.plotting.plot(detector_name, mantid.mtd[workspace])
        # for element, colour in zip(["Zn", "Cu"], ["b", "r"]):
        for element, colour in zip(["Zn"], ["b"]):
            data = self.ptable.element_data(element)["Primary"]
            for peak_type in self.ptable.element_data(element)["Primary"]:
                self.plotting.add_vline(
                    detector_name, data[peak_type], 0, 1, color=colour)
                #subplot.text(data[peak_type], 0.95*_max, element, fontsize=12, ha="center")

    def add_plot(self):
        name = "Plot {}".format(time())
        subplot = self.plotting.add_subplot(name)
        self.plotting.call_plot_method(name, subplot.set_title, name)
        plot1 = mantid.CreateSampleWorkspace(OutputWorkspace=str(time()))
        self.plotting.plot(name, plot1)
        plot2 = mantid.Plus(plot1, plot1, OutputWorkspace=str(time()))
        self.plotting.plot(name, plot2)
        self.plotting.add_hline(name, 0.06, 0, 1)
        self.plotting.add_vline(name, 10100, 0, 1)

    def del_plot(self):
        to_del = random.choice(self.plotting.get_subplots().keys())
        self.plotting.remove_subplot(to_del)

    def spinbox_changed(self, val):
        print("SpinBox Value Changed: {}".format(val))

    def spinbox_submit(self, val):
        print("SpinBox Submitted: {}".format(val))


def qapp():
    if QtGui.QApplication.instance():
        _app = QtGui.QApplication.instance()
    else:
        _app = QtGui.QApplication(sys.argv)
    return _app


app = qapp()
try:
    window = ElementalAnalysisGui()
    window.show()
    app.exec_()
except RuntimeError as error:
    message_box.warning(str(error))
