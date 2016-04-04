import unittest

import mantid
from mantid.api import (AlgorithmManager, mtd)
from SANSDataFileInfo import SANSDataFileInfo


SANS_DATA_FILE_INFO_SAMPLE_DATA_SET_NAME = "test_workspace_sans_data_file_info"

def load_sample_data_set():
    alg = AlgorithmManager.createUnmanaged("Load")
    alg.initialize()
    alg.setChild(True)
    alg.setProperty("Filename", "LOQ48127.nxs")
    alg.setProperty("OutputWorkspace", SANS_DATA_FILE_INFO_SAMPLE_DATA_SET_NAME)
    alg.execute()
    ws = alg.getProperty("OutputWorkspace").value
    mtd.addOrReplace(SANS_DATA_FILE_INFO_SAMPLE_DATA_SET_NAME, ws)
    return SANS_DATA_FILE_INFO_SAMPLE_DATA_SET_NAME

def delete_sample_data_set():
    if mtd.doesExist(SANS_DATA_FILE_INFO_SAMPLE_DATA_SET_NAME):
        mtd.remove(SANS_DATA_FILE_INFO_SAMPLE_DATA_SET_NAME)


class TestSANSDataFileInfo(unittest.TestCase):
    def test_that_data_which_is_loaded_is_not_reloaded(self):
        # Arrange
        ws_name = load_sample_data_set()
        ws_info = SANSDataFileInfo(ws_name)
        # Act + Assert
        self.assertTrue(ws_info.is_already_loaded(ws_name), "Should detect as being already loaded")
        # Clean up
        delete_sample_data_set()

    def test_that_data_which_is_not_loaded_is_loaded(self):
        # Arrange
        ws_name = load_sample_data_set()
        ws_info = SANSDataFileInfo(ws_name)
        ws_dummy_name = "test_name"
        # Act + Assert
        self.assertFalse(ws_info.is_already_loaded(ws_dummy_name), "Should detect has not been loaded")
        # Clean up
        delete_sample_data_set()

    def test_that_data_which_has_changed_its_size_is_reloaded(self):
        # Arrange
        import time
        time.sleep(10)
        ws_name = load_sample_data_set()
        ws_info = SANSDataFileInfo(ws_name)
        ws = mtd[ws_name]
        ws = ws + ws
        # Act + Assert
        self.assertFalse(ws_info.is_already_loaded(ws_name), "Should detect changed data set has to be reloaded")
        # Clean up
        delete_sample_data_set()

if __name__ == "__main__":
    unittest.main()
