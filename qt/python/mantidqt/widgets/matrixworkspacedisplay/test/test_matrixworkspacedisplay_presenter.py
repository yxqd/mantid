# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantid workbench.
#
#
from __future__ import (absolute_import, division, print_function)

import unittest

from mock import Mock

from mantidqt.widgets.matrixworkspacedisplay.presenter import MatrixWorkspaceDisplay
from mantidqt.widgets.matrixworkspacedisplay.test_helpers.matrixworkspacedisplay_common import MockQModelIndex, \
    MockWorkspace
from mantidqt.widgets.matrixworkspacedisplay.test_helpers.mock_matrixworkspacedisplay import \
    MockMatrixWorkspaceDisplayView, MockQTableView


class MatrixWorkspaceDisplayPresenterTest(unittest.TestCase):
    def assertNotCalled(self, mock):
        self.assertEqual(0, mock.call_count)

    def test_setup_table(self):
        ws = MockWorkspace()
        view = MockMatrixWorkspaceDisplayView()
        MatrixWorkspaceDisplay(ws, view=view)
        self.assertEqual(3, view.set_context_menu_actions.call_count)
        self.assertEqual(1, view.set_model.call_count)

    def test_action_copy_spectrum_values(self):
        ws = MockWorkspace()
        view = MockMatrixWorkspaceDisplayView()
        presenter = MatrixWorkspaceDisplay(ws, view=view)

        mock_table = MockQTableView()

        # two rows are selected in different positions
        mock_indexes = [MockQModelIndex(0, 1), MockQModelIndex(3, 1)]
        mock_table.mock_selection_model.selectedRows = Mock(return_value=mock_indexes)

        mock_read = Mock(return_value=[43, 99])
        presenter._get_ws_read_from_type = Mock(return_value=mock_read)
        expected_string = "43\t99\n43\t99"

        presenter.action_copy_spectrum_values(mock_table)

        mock_table.selectionModel.assert_called_once_with()
        mock_table.mock_selection_model.hasSelection.assert_called_once_with()
        view.copy_to_clipboard.assert_called_once_with(expected_string)
        view.show_mouse_toast.assert_called_once_with(MatrixWorkspaceDisplay.COPY_SUCCESSFUL_MESSAGE)

    def test_action_copy_spectrum_values_no_selection(self):
        ws = MockWorkspace()
        view = MockMatrixWorkspaceDisplayView()
        presenter = MatrixWorkspaceDisplay(ws, view=view)

        mock_table = MockQTableView()
        mock_table.mock_selection_model.hasSelection = Mock(return_value=False)
        mock_table.mock_selection_model.selectedRows = Mock()

        presenter.action_copy_spectrum_values(mock_table)

        mock_table.selectionModel.assert_called_once_with()
        mock_table.mock_selection_model.hasSelection.assert_called_once_with()
        # the action should never look for rows if there is no selection
        self.assertNotCalled(mock_table.mock_selection_model.selectedRows)
        view.show_mouse_toast.assert_called_once_with(MatrixWorkspaceDisplay.NO_SELECTION_MESSAGE)

    def test_action_copy_bin_values(self):
        ws = MockWorkspace()
        view = MockMatrixWorkspaceDisplayView()
        presenter = MatrixWorkspaceDisplay(ws, view=view)
        mock_table = MockQTableView()

        # two columns are selected at different positions
        mock_indexes = [MockQModelIndex(0, 0), MockQModelIndex(0, 3)]
        mock_table.mock_selection_model.selectedColumns = Mock(return_value=mock_indexes)
        # change the mock ws to have 3 histograms
        ws.getNumberHistograms = Mock(return_value=3)

        mock_read = Mock(return_value=[83, 11, 33, 70])
        presenter._get_ws_read_from_type = Mock(return_value=mock_read)
        expected_string = "83\t70\n83\t70\n83\t70"

        presenter.action_copy_bin_values(mock_table)

        mock_table.selectionModel.assert_called_once_with()
        view.copy_to_clipboard.assert_called_once_with(expected_string)
        mock_table.mock_selection_model.hasSelection.assert_called_once_with()
        view.show_mouse_toast.assert_called_once_with(MatrixWorkspaceDisplay.COPY_SUCCESSFUL_MESSAGE)

    def test_action_copy_bin_values_no_selection(self):
        ws = MockWorkspace()
        view = MockMatrixWorkspaceDisplayView()
        presenter = MatrixWorkspaceDisplay(ws, view=view)

        mock_table = MockQTableView()
        mock_table.mock_selection_model.hasSelection = Mock(return_value=False)
        mock_table.mock_selection_model.selectedColumns = Mock()

        presenter.action_copy_bin_values(mock_table)

        mock_table.selectionModel.assert_called_once_with()
        mock_table.mock_selection_model.hasSelection.assert_called_once_with()
        # the action should never look for rows if there is no selection
        self.assertNotCalled(mock_table.mock_selection_model.selectedColumns)
        view.show_mouse_toast.assert_called_once_with(MatrixWorkspaceDisplay.NO_SELECTION_MESSAGE)

    def test_action_copy_cell(self):
        ws = MockWorkspace()
        view = MockMatrixWorkspaceDisplayView()
        presenter = MatrixWorkspaceDisplay(ws, view=view)
        mock_table = MockQTableView()

        # two columns are selected at different positions
        mock_index = MockQModelIndex(None, None)
        mock_table.mock_selection_model.currentIndex = Mock(return_value=mock_index)

        presenter.action_copy_cells(mock_table)

        mock_table.selectionModel.assert_called_once_with()
        self.assertEqual(1, view.copy_to_clipboard.call_count)
        self.assertEqual(9, mock_index.sibling.call_count)
        view.show_mouse_toast.assert_called_once_with(MatrixWorkspaceDisplay.COPY_SUCCESSFUL_MESSAGE)

    def test_action_copy_cell_no_selection(self):
        ws = MockWorkspace()
        view = MockMatrixWorkspaceDisplayView()
        presenter = MatrixWorkspaceDisplay(ws, view=view)
        mock_table = MockQTableView()
        mock_table.mock_selection_model.hasSelection = Mock(return_value=False)

        presenter.action_copy_cells(mock_table)

        mock_table.selectionModel.assert_called_once_with()
        mock_table.mock_selection_model.hasSelection.assert_called_once_with()
        view.show_mouse_toast.assert_called_once_with(MatrixWorkspaceDisplay.NO_SELECTION_MESSAGE)

        self.assertNotCalled(view.copy_to_clipboard)

    def common_setup_action_plot(self, table_has_selection=True):
        mock_ws = MockWorkspace()
        mock_view = MockMatrixWorkspaceDisplayView()
        mock_plotter = Mock()
        presenter = MatrixWorkspaceDisplay(mock_ws, plot=mock_plotter, view=mock_view)

        # monkey-patch the spectrum plot label to count the number of calls
        presenter.model.get_spectrum_plot_label = Mock()
        presenter.model.get_bin_plot_label = Mock()

        mock_table = MockQTableView()
        # configure the mock return values
        mock_table.mock_selection_model.hasSelection = Mock(return_value=table_has_selection)
        return mock_plotter, mock_table, mock_view, presenter

    def setup_mock_selection(self, mock_table, num_selected_rows=None, num_selected_cols=None):
        """
        :type mock_table: MockQTableView
        :type num_selected_rows: int|None
        :type num_selected_cols: int|None
        """
        mock_selected = []
        if num_selected_rows is not None:
            for i in range(num_selected_rows):
                mock_selected.append(MockQModelIndex(i, 1))
            mock_table.mock_selection_model.selectedRows = Mock(return_value=mock_selected)
            mock_table.mock_selection_model.selectedColumns = Mock()
        elif num_selected_cols is not None:
            for i in range(num_selected_cols):
                mock_selected.append(MockQModelIndex(1, i))
            mock_table.mock_selection_model.selectedRows = Mock()
            mock_table.mock_selection_model.selectedColumns = Mock(return_value=mock_selected)
        else:
            mock_table.mock_selection_model.selectedRows = Mock()
            mock_table.mock_selection_model.selectedColumns = Mock()
        return mock_selected

    def test_action_plot_spectrum_plot_many_confirmed(self):
        mock_plot, mock_table, mock_view, presenter = self.common_setup_action_plot()
        num_selected_rows = MatrixWorkspaceDisplay.NUM_SELECTED_FOR_CONFIRMATION + 1

        self.setup_mock_selection(mock_table, num_selected_rows)

        # The a lot of things to plot message will show, set that the user will CONFIRM the plot
        # meaning the rest of the function will execute as normal
        mock_view.ask_confirmation = Mock(return_value=True)

        presenter.action_plot_spectrum(mock_table)

        mock_table.selectionModel.assert_called_once_with()
        mock_table.mock_selection_model.hasSelection.assert_called_once_with()
        mock_table.mock_selection_model.selectedRows.assert_called_once_with()

        mock_view.ask_confirmation.assert_called_once_with(
            MatrixWorkspaceDisplay.A_LOT_OF_THINGS_TO_PLOT_MESSAGE.format(num_selected_rows))

        self.assertNotCalled(mock_table.mock_selection_model.selectedColumns)
        self.assertEqual(1, mock_plot.call_count)

    def test_action_plot_spectrum_plot_many_denied(self):
        mock_plot, mock_table, mock_view, presenter = self.common_setup_action_plot()
        num_selected_rows = MatrixWorkspaceDisplay.NUM_SELECTED_FOR_CONFIRMATION + 1

        # return value unused as most of the function being tested is not executed
        self.setup_mock_selection(mock_table, num_selected_rows)

        # The a lot of things to plot message will show, set that the user will DENY the plot
        # meaning the rest of the function will NOT EXECUTE AT ALL
        mock_view.ask_confirmation = Mock(return_value=False)

        presenter.action_plot_spectrum(mock_table)

        mock_table.selectionModel.assert_called_once_with()
        mock_table.mock_selection_model.hasSelection.assert_called_once_with()
        mock_table.mock_selection_model.selectedRows.assert_called_once_with()

        mock_view.ask_confirmation.assert_called_once_with(
            MatrixWorkspaceDisplay.A_LOT_OF_THINGS_TO_PLOT_MESSAGE.format(num_selected_rows))

        self.assertNotCalled(mock_table.mock_selection_model.selectedColumns)
        self.assertNotCalled(mock_plot)

    def test_action_plot_spectrum_no_selection(self):
        mock_plot, mock_table, mock_view, presenter = self.common_setup_action_plot(table_has_selection=False)

        mock_table.mock_selection_model.selectedRows = Mock()
        mock_table.mock_selection_model.selectedColumns = Mock()

        presenter.action_plot_spectrum(mock_table)

        mock_view.show_mouse_toast.assert_called_once_with(MatrixWorkspaceDisplay.NO_SELECTION_MESSAGE)
        mock_table.selectionModel.assert_called_once_with()
        mock_table.mock_selection_model.hasSelection.assert_called_once_with()

        self.assertNotCalled(mock_table.mock_selection_model.selectedRows)
        self.assertNotCalled(mock_table.mock_selection_model.selectedColumns)
        self.assertNotCalled(mock_plot)

    def test_action_plot_bin_plot_many_confirmed(self):
        mock_plot, mock_table, mock_view, presenter = self.common_setup_action_plot()
        num_selected_cols = MatrixWorkspaceDisplay.NUM_SELECTED_FOR_CONFIRMATION + 1
        self.setup_mock_selection(mock_table, num_selected_cols=num_selected_cols)
        mock_view.ask_confirmation = Mock(return_value=True)

        presenter.action_plot_bin(mock_table)

        mock_table.selectionModel.assert_called_once_with()
        mock_table.mock_selection_model.hasSelection.assert_called_once_with()
        mock_table.mock_selection_model.selectedColumns.assert_called_once_with()
        self.assertNotCalled(mock_table.mock_selection_model.selectedRows)

        mock_view.ask_confirmation.assert_called_once_with(
            MatrixWorkspaceDisplay.A_LOT_OF_THINGS_TO_PLOT_MESSAGE.format(num_selected_cols))
        self.assertEqual(1, mock_plot.call_count)

    def test_action_plot_bin_plot_many_denied(self):
        mock_plot, mock_table, mock_view, presenter = self.common_setup_action_plot()
        num_selected_cols = MatrixWorkspaceDisplay.NUM_SELECTED_FOR_CONFIRMATION + 1

        # return value unused as most of the function being tested is not executed
        self.setup_mock_selection(mock_table, num_selected_cols=num_selected_cols)

        mock_view.ask_confirmation = Mock(return_value=False)

        presenter.action_plot_bin(mock_table)

        mock_table.selectionModel.assert_called_once_with()
        mock_table.mock_selection_model.hasSelection.assert_called_once_with()
        mock_table.mock_selection_model.selectedColumns.assert_called_once_with()

        mock_view.ask_confirmation.assert_called_once_with(
            MatrixWorkspaceDisplay.A_LOT_OF_THINGS_TO_PLOT_MESSAGE.format(num_selected_cols))

        self.assertNotCalled(mock_table.mock_selection_model.selectedRows)
        self.assertNotCalled(mock_plot)

    def test_action_plot_bin_no_selection(self):
        mock_plot, mock_table, mock_view, presenter = self.common_setup_action_plot(table_has_selection=False)
        self.setup_mock_selection(mock_table, num_selected_rows=None, num_selected_cols=None)

        presenter.action_plot_bin(mock_table)

        mock_view.show_mouse_toast.assert_called_once_with(MatrixWorkspaceDisplay.NO_SELECTION_MESSAGE)
        mock_table.selectionModel.assert_called_once_with()
        mock_table.mock_selection_model.hasSelection.assert_called_once_with()

        self.assertNotCalled(mock_table.mock_selection_model.selectedRows)
        self.assertNotCalled(mock_table.mock_selection_model.selectedColumns)
        self.assertNotCalled(mock_plot)


if __name__ == '__main__':
    unittest.main()
