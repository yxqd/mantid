from __future__ import (absolute_import, division, print_function)

from Muon import transform_presenter
from Muon import model_constructor

from Muon import transform_view

from Muon import help_view
from Muon import help_presenter


###### need to make a tab view constructora -> a dict to make it easier


class TabPresenter(object):
    """
    The widget for controlling which method to display
    in the transformation tab
    """
    def __init__(self,view,tabViews,groupedModels):
        self.view=view

        # create presenters for all of the tabs
        self.transformPresenter =transform_presenter.TransformPresenter(tabViews,groupedModels)
        helpView = help_view.HelpView()
        self.helpPresenter =help_presenter.HelpPresenter(helpView)
        #self.transformPresenter =transform_presenter.TransformPresenter(tabView.getTab("transform"),groupedModels)

    
        # add each widget to taab
        #self.view.addTab(tabView.getTab("transform"),"transform")
        self.view.addTab(tabViews,"transform")
        self.view.addTab(helpView,"Help")

      
        # make into tabs
        self.view.makeTabs()


    def close(self):
        self.transformPresenter.close()

