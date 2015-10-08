

class ProcessingViewPresenter(object):

    def __init__(self, processing_view):
        self._processing_view = processing_view

    def show(self):
        self._processing_view.show_empty()
