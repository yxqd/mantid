
import unittest

can_mock = True
try:
    import mock
except ImportError:
    can_mock = False

from main_view_presenter import MainViewPresenter
from view.main_view import MainView
from command import Command

class MainViewPresenterTest(unittest.TestCase):

    def test_construction(self):

        if not can_mock:
            return

        main_view = mock.create_autospec(MainView)

        presenter = MainViewPresenter(main_view)

    def test_process_all(self):
        if not can_mock:
            return

        main_view = mock.create_autospec(MainView)

        presenter = MainViewPresenter(main_view)

        presenter.notify(Command.ProcessAll)

        self.assertEqual(1, main_view.get_run_table.call_count)
        self.assertEqual(1, main_view.set_processing.call_count)



if __name__ == '__main__':
    unittest.main()
