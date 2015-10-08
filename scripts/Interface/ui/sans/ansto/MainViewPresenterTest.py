
import unittest

can_mock = True
try:
    import mock
except ImportError:
    can_mock = False

from presenter.main_view_presenter import MainViewPresenter
from view.main_view import MainView

class MainViewPresenterTest(unittest.TestCase):

    def test_construction(self):
        if not can_mock:
            return

        main_view = mock.create_autospec(MainView)

        presenter = MainViewPresenter(main_view)

        self.assertEqual(1, main_view.get_processing_view.call_count, msg="Should be fetching the processing view first-off")
        self.assertEqual(1, main_view.select_processing_view.call_count, msg="Should be selecting the processing view")


if __name__ == '__main__':
    unittest.main()
