from __future__ import (absolute_import, division, print_function)

from Muon import transform_presenter
from Muon import help_presenter

# only in here as a dummy tab
from Muon import help_view


class TabPresenter(object):
    """
    The widget for controlling which method to display
    in the transformation tab
    """
    def __init__(self,view,tabViews,groupedModels):
        self.view=view
        self.createPresenters(tabViews,groupedModels)
      
        # make into tabs
        self.view.makeTabs()


    """
    Only the following code would need updating for a new GUI
    """
    def createPresenters(self,tabViews,groupedModels):
        """
        This code may need updating for a different GUI
        It constructs the different presenters that are contained 
        within the main window
        """
        # create presenters for all of the tabs
        self.transformPresenter =transform_presenter.TransformPresenter(tabViews,groupedModels)
        helpView = help_view.HelpView()
        self.helpPresenter =help_presenter.HelpPresenter(helpView)

        # add each widget to tab
        self.view.addTab(tabViews,"transform")
        self.view.addTab(helpView,"Help")


    def close(self):
        """
        This code may need updating for a different GUI
        It contains the close command for all of the presenters
        contained in the main window
        """
        self.transformPresenter.close()
        #self.helpPresenter.close()
