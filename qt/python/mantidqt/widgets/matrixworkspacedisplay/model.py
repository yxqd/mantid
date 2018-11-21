# coding=utf-8
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

from mantid.api import MatrixWorkspace
from mantid.dataobjects import EventWorkspace, Workspace2D
from mantidqt.widgets.matrixworkspacedisplay.table_view_model import MatrixWorkspaceTableViewModel, \
    MatrixWorkspaceTableViewModelType
from mantidqt.widgets.matrixworkspacedisplay.test_helpers.matrixworkspacedisplay_common import MockWorkspace


class MatrixWorkspaceDisplayModel(object):
    SPECTRUM_PLOT_LEGEND_STRING = '{}-{}'
    BIN_PLOT_LEGEND_STRING = '{}-bin-{}'

    def __init__(self, ws):
        if not isinstance(ws, MatrixWorkspace) \
                and not isinstance(ws, Workspace2D) \
                and not isinstance(ws, EventWorkspace) \
                and not isinstance(ws, MockWorkspace):
            raise ValueError("The workspace type is not supported: {0}".format(type(ws)))

        self._ws = ws

    def get_name(self):
        return self._ws.name()

    def get_item_model(self):
        return MatrixWorkspaceTableViewModel(self._ws, MatrixWorkspaceTableViewModelType.x), \
               MatrixWorkspaceTableViewModel(self._ws, MatrixWorkspaceTableViewModelType.y), \
               MatrixWorkspaceTableViewModel(self._ws, MatrixWorkspaceTableViewModelType.e)
