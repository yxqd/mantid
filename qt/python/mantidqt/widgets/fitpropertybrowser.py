# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#
#
from __future__ import (absolute_import, unicode_literals)

from qtpy.QtCore import Signal

from mantidqt.utils.qt import import_qt


BaseBrowser = import_qt('.._common', 'mantidqt.widgets', 'FitPropertyBrowser')


class FitPropertyBrowser(BaseBrowser):

    closing = Signal()

    def __init__(self, parent=None):
        super(FitPropertyBrowser, self).__init__(parent)
        self.init()

    def closeEvent(self, event):
        self.closing.emit()
        BaseBrowser.closeEvent(self, event)
