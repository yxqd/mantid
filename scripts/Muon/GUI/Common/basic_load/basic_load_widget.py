from __future__ import (absolute_import, division, print_function)

from Muon.GUI.Common.basic_load.basic_load_view import BasicLoadView
from Muon.GUI.Common.basic_load.basic_load_presenter import BasicLoadPresenter


class BasicLoadWidget(object):

    def __init__(self, name, parent=None):
        view = BasicLoadView(name, parent)
        model = None
        self.presenter = BasicLoadPresenter(view, model)

    def getPresenter(self):
        return self.presenter

    @property
    def widget(self):
        return self.presenter.widget

