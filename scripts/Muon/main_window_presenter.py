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
        self.groupedModels=groupedModels
        #view =transform_view.TransformView(self.groupedViews,parent)
        self.Transformpresenter =transform_presenter.TransformPresenter(self.view.transformView,self.groupedModels)

        groupedModels= model_constructor.ModelConstructor(True)
        self.tabPresenter = tab_presenter.TabPresenter(self.view.tabView,self.view.transformView,groupedModels)
