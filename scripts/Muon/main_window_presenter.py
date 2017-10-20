from __future__ import (absolute_import, division, print_function)


from Muon import tab_presenter


class MainWindowPresenter(object):
    """
    This is the presenter for the main GUI.
    It connstructs the oberall look of the GUI
    """
    def __init__(self,view,groupedModels):
        self.view=view
        self.createPresenters(groupedModels)


    """
    Only the following code would need updating for a new GUI
    """
    def createPresenters(self,groupedModels):
        """
        This code may need updating for a different GUI
        It constructs the different presenters that are contained 
        within the main window
        """
        childViews = self.view.getViews()
 
        self.tabPresenter = tab_presenter.TabPresenter(childViews["Tab"],childViews["Tab"].transformView,groupedModels)

    def close(self):
        """
        This code may need updating for a different GUI
        It contains the close command for all of the presenters
        contained in the main window
        """
        self.tabPresenter.close()
