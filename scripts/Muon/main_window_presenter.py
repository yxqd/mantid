from __future__ import (absolute_import, division, print_function)


from Muon import model_constructor
from Muon import transform_presenter
from Muon import tab_presenter
from Muon import transform_view


class MainWindowPresenter(object):
    """
    This is the presenter for the main GUI.
    It connstructs the oberall look of the GUI
    """
    def __init__(self,view,groupedModels):
        self.view=view

        groupedModels= model_constructor.ModelConstructor(True)
        self.tabPresenter = tab_presenter.TabPresenter(self.view.tabView,self.view.tabView.transformView,groupedModels)
