from __future__ import (absolute_import, division, print_function)


from Muon import model_constructor
from Muon import transform_presenter
from Muon import transform_view


class MainWindowPresenter(object):
    """
    This is the presenter for the main GUI.
    It connstructs the oberall look of the GUI
    """
    def __init__(self,view,groupedViews,groupedModels):
        self.view=view
        self.groupedViews=groupedViews
        self.groupedModels=groupedModels
        #view =transform_view.TransformView(self.groupedViews,parent)
        self.Transformpresenter =transform_presenter.TransformPresenter(self.view,self.groupedModels)

 
