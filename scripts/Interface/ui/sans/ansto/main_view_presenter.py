from command import Command

class MainViewPresenter(object):

    def __init__(self, main_view):
        self._main_view = main_view
        # Default view is the processing view

    def notify(self, command):
        if command == Command.ProcessAll:

            """
            TODO. Run the data processing algorithm over all runs

            1. Extract settings from settings tab (via settings_view_presenter TODO)
            2. For loop over all table inputs (runs)
            3. Load necessary files via Load()
            4. Run the DataProcessorAlgorithm
            5. Run any other algorithm step
            """
            table_of_runs = self._main_view.get_run_table()
            self._main_view.set_processing()
        else:
            # TODO more command driven processing options
            pass

