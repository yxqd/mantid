from __future__ import (absolute_import, division, print_function)


class DummyPresenter(object):

    def __init__(self,view,model):
        self._view = view
        self._model = model

    @property
    def widget(self):
        return self._view
