# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
from __future__ import (absolute_import, division, print_function)

from qtpy import QtWidgets, QtCore
from Muon.GUI.Common.utilities import table_utils
from Muon.GUI.Common.message_box import warning
from mantidqt.utils.qt import load_ui
from mantidqt.widgets.functionbrowser import FunctionBrowser

ui_fitting_tab, _ = load_ui(__file__, "fitting_tab.ui")


class FittingTabView(QtWidgets.QWidget, ui_fitting_tab):
    def __init__(self, parent=None):
        super(FittingTabView, self).__init__(parent)
        self.setupUi(self)
        self.setup_fit_options_table()

        self.function_browser = FunctionBrowser(self, True)
        self.function_browser_layout.addWidget(self.function_browser)

        self.increment_parameter_display_button.clicked.connect(self.increment_display_combo_box)
        self.decrement_parameter_display_button.clicked.connect(self.decrement_display_combo_box)

    def update_displayed_data_combo_box(self, data_list):
        name = self.parameter_display_combo.currentText()
        self.parameter_display_combo.clear()

        self.parameter_display_combo.addItems(data_list)

        index = self.parameter_display_combo.findText(name)

        if index != -1:
            self.parameter_display_combo.setCurrentIndex(index)
        else:
            self.parameter_display_combo.setCurrentIndex(0)

    def increment_display_combo_box(self):
        index = self.parameter_display_combo.currentIndex()
        count = self.parameter_display_combo.count()

        if index < count - 1:
            self.parameter_display_combo.setCurrentIndex(index + 1)
        else:
            self.parameter_display_combo.setCurrentIndex(0)

    def decrement_display_combo_box(self):
        index = self.parameter_display_combo.currentIndex()
        count = self.parameter_display_combo.count()

        if index != 0:
            self.parameter_display_combo.setCurrentIndex(index - 1)
        else:
            self.parameter_display_combo.setCurrentIndex(count - 1)

    def set_datasets_in_function_browser(self, data_set_name_list):
        number_of_data_sets = self.function_browser.getNumberOfDatasets()
        index_list = range(number_of_data_sets)
        self.function_browser.removeDatasets(index_list)

        self.function_browser.addDatasets(data_set_name_list)

    def update_with_fit_outputs(self, fit_function, output_status, output_chi_squared):
        self.function_browser.setFunction(str(fit_function))
        if output_status == 'success':
            self.fit_status_success_failure.setText('Success')
            self.fit_status_success_failure.setStyleSheet('color: green')
        else:
            self.fit_status_success_failure.setText('Failure: {}'.format(output_status))
            self.fit_status_success_failure.setStyleSheet('color: red')

        self.fit_status_chi_squared.setText('Chi squared: {}'.format(output_chi_squared))

    def update_global_fit_state(self, number_succesfully_fitted, total_number_fitted):
        self.global_fit_status_label.setText('Fit successful for {} of {} workspaces'.
                                             format(number_succesfully_fitted, total_number_fitted))

        if number_succesfully_fitted == total_number_fitted:
            self.global_fit_status_label.setStyleSheet('color: green')
        else:
            self.global_fit_status_label.setStyleSheet('color: red')

    def set_slot_for_select_workspaces_to_fit(self, slot):
        self.select_workspaces_to_fit_button.clicked.connect(slot)

    def set_slot_for_display_workspace_changed(self, slot):
        self.parameter_display_combo.currentIndexChanged.connect(slot)

    def set_slot_for_use_raw_changed(self, slot):
        self.fit_to_raw_data_checkbox.stateChanged.connect(slot)

    def set_slot_for_fit_type_changed(self, slot):
        self.single_fit_radio.toggled.connect(slot)
        self.simul_fit_radio.toggled.connect(slot)
        self.sequential_fit_radio.toggled.connect(slot)

    def set_slot_for_fit_button_clicked(self, slot):
        self.fit_button.clicked.connect(slot)

    def set_slot_for_start_x_updated(self, slot):
        self.time_start.editingFinished.connect(slot)

    def set_slot_for_end_x_updated(self, slot):
        self.time_end.editingFinished.connect(slot)

    @property
    def display_workspace(self):
        return str(self.parameter_display_combo.currentText())

    @property
    def fit_string(self):
        return str(self.function_browser.getFitFunctionString())

    @property
    def minimizer(self):
        return str(self.minimizer_combo.currentText())

    @property
    def start_time(self):
        return float(self.time_start.text())

    @start_time.setter
    def start_time(self, value):
        self.time_start.setText(str(value))

    @property
    def end_time(self):
        return float(self.time_end.text())

    @end_time.setter
    def end_time(self, value):
        self.time_end.setText(str(value))

    @property
    def evaluation_type(self):
        return str(self.evaluation_combo.currentText())

    @property
    def fit_type(self):
        if self.single_fit_radio.isChecked():
            return self.single_fit_radio.text()
        if self.simul_fit_radio.isChecked():
            return self.simul_fit_radio.text()
        if self.sequential_fit_radio.isChecked():
            return self.sequential_fit_radio.text()

    @property
    def fit_to_raw(self):
        return self.fit_to_raw_data_checkbox.isChecked()

    @fit_to_raw.setter
    def fit_to_raw(self, value):
        state = QtCore.Qt.Checked if value else QtCore.Qt.Unchecked
        self.fit_to_raw_data_checkbox.setCheckState(state)

    def warning_popup(self, message):
        warning(message, parent=self)

    def get_index_for_start_end_times(self):
        if self.fit_type == 'Single Fit':
            return 0

        current_index = self.parameter_display_combo.currentIndex()
        current_index = current_index if current_index != -1 else 0

        return current_index

    def setup_fit_options_table(self):
        self.fit_options_table.setRowCount(5)
        self.fit_options_table.setColumnCount(2)
        self.fit_options_table.setColumnWidth(0, 300)
        self.fit_options_table.setColumnWidth(1, 300)
        self.fit_options_table.verticalHeader().setVisible(False)
        self.fit_options_table.horizontalHeader().setStretchLastSection(True)
        self.fit_options_table.setHorizontalHeaderLabels(
            ("Property;Value").split(";"))

        table_utils.setRowName(self.fit_options_table, 0, "Time Start")
        self.time_start = table_utils.addDoubleToTable(self.fit_options_table, 0.0, 0, 1)

        table_utils.setRowName(self.fit_options_table, 1, "Time End")
        self.time_end = table_utils.addDoubleToTable(self.fit_options_table, 15.0, 1, 1)

        table_utils.setRowName(self.fit_options_table, 2, "Minimizer")
        self.minimizer_combo = table_utils.addComboToTable(self.fit_options_table, 2, [])

        self.minimizer_combo.addItems(
            ['Levenberg-Marquardt', 'BFGS', 'Conjugate gradient (Fletcher-Reeves imp.)', 'Conjugate gradient (Polak-Ribiere imp.)',
             'Damped GaussNewton', 'FABADA', 'Levenberg-MarquardtMD', 'Simplex',
             'SteepestDescent', 'Trust Region'])

        table_utils.setRowName(self.fit_options_table, 3, "Fit To Raw Data")
        self.fit_to_raw_data_checkbox = table_utils.addCheckBoxWidgetToTable(
            self.fit_options_table, True, 3)

        table_utils.setRowName(self.fit_options_table, 4, "Evaluate Function As")
        self.evaluation_combo = table_utils.addComboToTable(self.fit_options_table, 4, ['CentrePoint', 'Histogram'])
