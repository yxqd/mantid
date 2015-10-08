import unittest

can_mock = True
try:
    import mock
except ImportError:
    can_mock = False

from presenter.processing_view_presenter import ProcessingViewPresenter
from view.processing_view import ProcessingView

class ProcessingViewPresenterTest(unittest.TestCase):

    def test_construction(self):
        if not can_mock:
            return

        processing_view = mock.create_autospec(ProcessingView)

        presenter = ProcessingViewPresenter(processing_view)


    def test_show(self):
        if not can_mock:
            return

        processing_view = mock.create_autospec(ProcessingView)

        presenter = ProcessingViewPresenter(processing_view)
        presenter.show()

        self.assertEqual(1, processing_view.show_empty.call_count, msg="Should be instructed to show nothing initially")

if __name__ == '__main__':
    unittest.main()