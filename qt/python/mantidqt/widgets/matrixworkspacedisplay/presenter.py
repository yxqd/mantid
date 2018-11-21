# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantid workbench.
#
#
from __future__ import absolute_import, division, print_function

from mantid.plots.utility import MantidAxType
from mantidqt.widgets.matrixworkspacedisplay.table_view_model import MatrixWorkspaceTableViewModelType
from .model import MatrixWorkspaceDisplayModel
from .view import MatrixWorkspaceDisplayView


class MatrixWorkspaceDisplay(object):
    NO_SELECTION_MESSAGE = "No selection"
    COPY_SUCCESSFUL_MESSAGE = "Copy Successful"
    A_LOT_OF_THINGS_TO_PLOT_MESSAGE = "You selected {} spectra to plot. Are you sure you want to plot that many?"
    NUM_SELECTED_FOR_CONFIRMATION = 10

    def __init__(self, ws, plot=None, parent=None, model=None, view=None):
        # Create model and view, or accept mocked versions
        self.model = model if model else MatrixWorkspaceDisplayModel(ws)
        self.view = view if view else MatrixWorkspaceDisplayView(self,
                                                                 parent,
                                                                 self.model.get_name())

        self.plot = plot
        self.setup_tables()

        self.view.set_context_menu_actions(self.view.table_y)
        self.view.set_context_menu_actions(self.view.table_x)
        self.view.set_context_menu_actions(self.view.table_e)

    def setup_tables(self):
        # unpacks the list of models returned from getItemModel
        self.view.set_model(*self.model.get_item_model())

    def action_copy_spectrum_values(self, table):
        """
        Copies the values selected by the user to the system's clipboard

        :param table: Table from which the selection will be read
        :param ws_read: The workspace read function, that is used to access the data directly
        """
        selection_model = table.selectionModel()
        if not selection_model.hasSelection():
            self.show_no_selection_to_copy_toast()
            return
        ws_read = self._get_ws_read_from_type(table.model().type)
        selected_rows = selection_model.selectedRows()  # type: list
        row_data = []

        for index in selected_rows:
            row = index.row()
            data = "\t".join(map(str, ws_read(row)))

            row_data.append(data)

        self.view.copy_to_clipboard("\n".join(row_data))
        self.show_successful_copy_toast()

    def show_no_selection_to_copy_toast(self):
        self.view.show_mouse_toast(self.NO_SELECTION_MESSAGE)

    def show_successful_copy_toast(self):
        self.view.show_mouse_toast(self.COPY_SUCCESSFUL_MESSAGE)

    def action_copy_bin_values(self, table):
        selection_model = table.selectionModel()
        if not selection_model.hasSelection():
            self.show_no_selection_to_copy_toast()
            return
        ws_read = self._get_ws_read_from_type(table.model().type)
        selected_columns = selection_model.selectedColumns()  # type: list

        # Qt gives back a QModelIndex, we need to extract the column from it
        num_rows = self.model._ws.getNumberHistograms()
        column_data = []
        for index in selected_columns:
            column = index.column()
            data = [str(ws_read(row)[column]) for row in range(num_rows)]
            column_data.append(data)

        all_string_rows = []
        for i in range(num_rows):
            # Appends ONE value from each COLUMN, this is because the final string is being built vertically
            # the noqa disables a 'data' variable redefined warning
            all_string_rows.append("\t".join([data[i] for data in column_data]))  # noqa: F812

        # Finally all rows are joined together with a new line at the end of each row
        final_string = "\n".join(all_string_rows)
        self.view.copy_to_clipboard(final_string)
        self.show_successful_copy_toast()

    def action_copy_cells(self, table):
        """
        :type table: QTableView
        :param table: The table from which the data will be copied.
        :return:
        """
        selectionModel = table.selectionModel()
        if not selectionModel.hasSelection():
            self.show_no_selection_to_copy_toast()
            return

        selection = selectionModel.selection()
        selectionRange = selection.first()

        top = selectionRange.top()
        bottom = selectionRange.bottom()
        left = selectionRange.left()
        right = selectionRange.right()

        data = []
        index = selectionModel.currentIndex()
        for i in range(top, bottom + 1):
            for j in range(left, right):
                data.append(index.sibling(i, j).data())
                data.append("\t")
            data.append(index.sibling(i, right).data())
            data.append("\n")

        # strip the string to remove the trailing new line
        self.view.copy_to_clipboard("".join(data).strip())
        self.show_successful_copy_toast()

    def _do_action_plot(self, table, axis, get_index, plot_errors=False):
        if self.plot is None:
            raise ValueError("Trying to do a plot, but no plotting class dependency was injected in the constructor")
        selection_model = table.selectionModel()
        if not selection_model.hasSelection():
            self.show_no_selection_to_copy_toast()
            return

        if axis == MantidAxType.SPECTRUM:
            selected = selection_model.selectedRows()  # type: list
        else:
            selected = selection_model.selectedColumns()  # type: list

        if len(selected) > self.NUM_SELECTED_FOR_CONFIRMATION and not self.view.ask_confirmation(
                self.A_LOT_OF_THINGS_TO_PLOT_MESSAGE.format(len(selected))):
            return

        plot_kwargs = {"capsize": 3} if plot_errors else {}
        plot_kwargs["axis"] = axis

        ws_list = [self.model._ws]
        self.plot(ws_list, wksp_indices=[get_index(index) for index in selected], errors=plot_errors,
                  plot_kwargs=plot_kwargs)

    def action_plot_spectrum(self, table):
        self._do_action_plot(table, MantidAxType.SPECTRUM, lambda index: index.row())

    def action_plot_spectrum_with_errors(self, table):
        self._do_action_plot(table, MantidAxType.SPECTRUM, lambda index: index.row(), plot_errors=True)

    def action_plot_bin(self, table):
        self._do_action_plot(table, MantidAxType.BIN, lambda index: index.column())

    def action_plot_bin_with_errors(self, table):
        self._do_action_plot(table, MantidAxType.BIN, lambda index: index.column(), plot_errors=True)

    def action_keypress_copy(self, table):
        selectionModel = table.selectionModel()
        if not selectionModel.hasSelection():
            self.show_no_selection_to_copy_toast()
            return

        if len(selectionModel.selectedRows()) > 0:
            self.action_copy_spectrum_values(table)
        elif len(selectionModel.selectedColumns()) > 0:
            self.action_copy_bin_values(table)
        else:
            self.action_copy_cells(table)

    def _get_ws_read_from_type(self, type):
        if type == MatrixWorkspaceTableViewModelType.y:
            return self.model._ws.readY
        elif type == MatrixWorkspaceTableViewModelType.x:
            return self.model._ws.readX
        elif type == MatrixWorkspaceTableViewModelType.e:
            return self.model._ws.readE
        else:
            raise ValueError("Unknown TableViewModel type {}".format(type))
