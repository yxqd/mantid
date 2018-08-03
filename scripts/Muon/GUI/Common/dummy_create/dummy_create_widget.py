from __future__ import (absolute_import, division, print_function)

from Muon.GUI.Common.dummy_create.dummy_create_view import DummyCreateView
from Muon.GUI.Common.dummy_create.dummy_create_presenter import DummyCreatePresenter
from Muon.GUI.Common.dummy_create.dummy_create_model import DummyCreateModel, DummyCreateWrapper

from PyQt4 import QtGui


class DummyCreateWidget(QtGui.QWidget):

    def __init__(self, parent=None,reporter=None):
        super(DummyCreateWidget, self).__init__(parent)
        view = DummyCreateView(parent)

        tmp = DummyCreateModel()
        model = DummyCreateWrapper(tmp)

        self._presenter = DummyCreatePresenter(view=view, alg=model,reporter=reporter)

    @property
    def presenter(self):
        return self._presenter

    @property
    def widget(self):
        return self._presenter.widget

