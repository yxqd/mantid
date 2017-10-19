from __future__ import (absolute_import, division, print_function)


class HelpPresenter(object):
    """
    This is the presenter for the help widget.
    It connects the view and model together and deals with
    logic.
    """
    def __init__(self,view):
        self.view=view
 
