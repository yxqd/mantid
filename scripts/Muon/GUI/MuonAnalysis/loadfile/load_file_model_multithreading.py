from __future__ import (absolute_import, division, print_function)

import os

import mantid.simpleapi as mantid
from mantid.kernel import ConfigService
from mantid import config as cf

from Muon.GUI.Common.muon_load_data import MuonLoadData
import Muon.GUI.Common.threading_manager as thread_manager


class BrowseFileWidgetModel(object):

    def __init__(self, loaded_data_store=MuonLoadData()):
        # Temporary list of filenames used for load thread
        self._filenames = []

        self._loaded_data_store = loaded_data_store

        self.thread_manager = None

    @property
    def loaded_filenames(self):
        return self._loaded_data_store.get_parameter("filename")

    @property
    def loaded_workspaces(self):
        return self._loaded_data_store.get_parameter("workspace")

    @property
    def loaded_runs(self):
        return self._loaded_data_store.get_parameter("run")

    # Used with load thread
    def output(self):
        pass

    # Used with load thread
    def cancel(self):
        pass

    # Used with load thread
    def loadData(self, filename_list):
        self._filenames = filename_list

    # Used with load thread
    def execute(self):
        failed_files = []
        for filename in self._filenames:
            try:
                ws, run = self.load_workspace_from_filename(filename)
            except ValueError:
                failed_files += [filename]
                continue
            self._loaded_data_store.add_data(run=run, workspace=ws, filename=filename)
        if failed_files:
            message = self.exception_message_for_failed_files(failed_files)
            raise ValueError(message)

    def add_thread_data(self):
        print("ADDING THREAD DATA : ")
        for res in self.thread_manager.results["results"]:
            self._loaded_data_store.add_data(run=res[2], workspace=res[1], filename=res[0])

    def load_workspace_from_filename(self, filename):
        print("Loading file : ", filename)
        try:
            workspace = mantid.Load(Filename=filename)
            run = int(workspace[0].getRunNumber())
        except ValueError as e:
            raise ValueError(e.args)
        return filename, workspace, run

    def load_with_multithreading(self, filenames, callback_finished, callback_progress, callback_exception):
        self.load_func = thread_manager.threading_decorator(self.load_workspace_from_filename)
        n_threads = min(2, len(filenames))
        self.thread_manager = thread_manager.WorkerManager(fn = self.load_workspace_from_filename, num_threads=n_threads,
                                                           callback_on_progress_update=callback_progress,
                                                           callback_on_thread_exception=callback_exception,
                                                           callback_on_threads_complete=callback_finished,
                                                           filename=filenames)
        print("STARTING MULTI THREADING")
        self.thread_manager.start()

    def exception_message_for_failed_files(self, failed_file_list):
        print("Exception in execute!")
        return "Could not load the following files : \n - " + "\n - ".join(failed_file_list)


    def clear(self):
        self._loaded_data_store.clear()

    def add_directories_to_config_service(self, file_list):
        print(file_list)
        dirs = [os.path.dirname(filename) for filename in file_list]
        dirs = [filename if os.path.isdir(filename) else "" for filename in dirs]
        dirs = list(set(dirs))
        if dirs:
            for dir in dirs:
                ConfigService.Instance().appendDataSearchDir(dir.encode('ascii', 'ignore'))
