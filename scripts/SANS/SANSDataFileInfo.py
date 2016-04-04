from mantid.simpleapi import *


class SANSDataFileInfo(object):
    def __init__(self, workspace_name):
        super(SANSDataFileInfo, self).__init__()
        self._workspace_name = workspace_name
        self._full_file_path = None
        self._history = None
        self._size = None

    def _initialize(self):
        if self._workspace_name is not None and mtd.doesExist(self._workspace_name):
            self._size = self.get_size_from_workspace(self._workspace_name)
        else:
            self._reset()

    def _reset(self):
        self._workspace_name = None
        self._full_file_path = None
        self._history = None
        self._size = None

    @staticmethod
    def get_size_from_workspace(workspace_name):
        if mtd.doesExist(workspace_name):
            ws = mtd.retrieve(workspace_name)
            return 0
            # TODO: Finish getting size
        else:
            return 0

    @staticmethod
    def get_history_from_workspace(workspace_name):
        if mtd.doesExist(workspace_name):
            ws = mtd.retrieve(workspace_name)
            return ""
            # TODO: Finish getting history
        else:
            return ""

    @staticmethod
    def get_full_file_path_from_workspace(workspace_name):
        if mtd.doesExist(workspace_name):
            ws = mtd.retrieve(workspace_name)
            return ""
            # TODO: Finish getting full file path
        else:
            return ""

    def is_already_loaded(self, requested_workspace_name):
        # If workspace name is None, or none is loaded or names don't match
        if (self._workspace_name is None or not mtd.doesExist(self._workspace_name) or
                not self.is_same_workspace_name(self._workspace_name, requested_workspace_name) or
                not self._workspace_is_still_the_same(requested_workspace_name)):
            return False
        else:
            return True

    def _workspace_is_still_the_same(self, requested_workspace_name):
        if self._workspace_name is None or not mtd.doesExist(self._workspace_name):
            return (self._are_same_file_paths(requested_workspace_name) and
                    self._are_same_compare_sizes(requested_workspace_name) and
                    self._are_same_compare_histories(requested_workspace_name))
        else:
            return True

    def _are_same_file_paths(self, requested_workspace_name):
        reference_file_path = self.get_full_file_path_from_workspace(requested_workspace_name)
        if reference_file_path == self._full_file_path:
            return True
        else:
            return False

    def _are_same_sizes(self, requested_workspace_name):
        reference_size = self.get_size_from_workspace(requested_workspace_name)
        if reference_size == self._size:
            return True
        else:
            return False

    def _are_same_histories(self, requested_workspace_name):
        reference_history = self.get_history_from_workspace(requested_workspace_name)
        if reference_history == self._history:
            return True
        else:
            return False

    @staticmethod
    def is_same_workspace_name(workspace_name_1, workspace_name_2):
        return workspace_name_1 == workspace_name_2

class SANSDataFileInfoStorage(object):
    def __init__(self):
        super(SANSDataFileInfoStorage, self).__init__()
        self.sample_sans = None
        self.sample_transmission = None
        self.sample_direct = None
        self.can_sans = None
        self.can_transmission = None
        self.can_direct = None

    def set_sample_sans(self, workspace_name):
        pass