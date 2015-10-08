from processing_view_presenter import ProcessingViewPresenter

class MainViewPresenter(object):

    def __init__(self, main_view):
        self._main_view = main_view
        # Default view is the processing view
        self._processing_presenter = ProcessingViewPresenter(main_view.get_processing_view())
        # Show the processing view
        self._main_view.select_processing_view()
        # Let the sub presenter show
        self._processing_presenter.show()
