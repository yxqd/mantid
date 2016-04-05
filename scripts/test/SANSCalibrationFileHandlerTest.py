import unittest

import mantid
from mantid.api import (AlgorithmManager, mtd)
from SANSCalibrationFileHandler import SANSCalibrationFileHandler
from SANSUtility import getFileAndName


# Files which are used to test the loading of the calibration files (this is all dummy)
def load_data_set(file_name, name):
    alg = AlgorithmManager.createUnmanaged("Load")
    alg.initialize()
    alg.setChild(True)
    alg.setProperty("Filename", file_name)
    alg.setProperty("OutputWorkspace", name)
    alg.execute()
    ws = alg.getProperty("OutputWorkspace").value
    mtd.addOrReplace(name, ws)

def delete_data_set(name):
    if mtd.doesExist(name):
        mtd.remove(name)


class TestCalibratitionFileHandler(unittest.TestCase):
    def test_that_initially_file_is_loaded(self):
        # Arrange
        file_name = "LOQ48127.nxs"
        calibration_file_path, calibration_file_name = getFileAndName(file_name)
        calibration_file_handler = SANSCalibrationFileHandler()

        # Act
        ws = calibration_file_handler.provide_calibration_file(calibration_file_name, calibration_file_path)

        # Assert
        self.assertTrue(ws is not None, "A workspace should have been loaded")

        # Clean up
        delete_data_set(calibration_file_name)

    def test_that_already_loaded_file_is_not_reloaded(self):
        # Arrange
        file_name = "LOQ48127.nxs"
        calibration_file_path, calibration_file_name = getFileAndName(file_name)
        calibration_file_handler = SANSCalibrationFileHandler()
        calibration_file_handler.provide_calibration_file(calibration_file_name, calibration_file_path)

        # Act + Assert
        self.assertTrue(calibration_file_handler.requires_no_reload(calibration_file_name, calibration_file_path),
                         "The data set should not be reloaded")

        # Clean up
        delete_data_set(calibration_file_name)

    def test_that_new_calibration_file_is_loaded_if_name_is_different(self):
        # Arrange
        file_name = "LOQ48127.nxs"
        calibration_file_path, calibration_file_name = getFileAndName(file_name)
        calibration_file_handler = SANSCalibrationFileHandler()
        ws = calibration_file_handler.provide_calibration_file(calibration_file_name, calibration_file_path)

        new_file_name = "LOQ48127np.nxs"
        new_calibration_file_path, new_calibration_file_name = getFileAndName(new_file_name)

        # Act + Assert
        self.assertFalse(calibration_file_handler.requires_no_reload(new_calibration_file_name, new_calibration_file_path),
                         "The data set should be reloaded")

        # Clean up
        delete_data_set(calibration_file_name)
        delete_data_set(new_calibration_file_name)

    def test_that_new_calibration_file_is_loaded_if_calibration_file_has_been_updated(self):
        # Arrange
        file_name = "LOQ48127.nxs"
        calibration_file_path, calibration_file_name = getFileAndName(file_name)
        calibration_file_handler = SANSCalibrationFileHandler()
        ws = calibration_file_handler.provide_calibration_file(calibration_file_name, calibration_file_path)

        # Act + Assert
        calibration_file_handler.calibration_was_updated = True
        self.assertFalse(calibration_file_handler.requires_no_reload(calibration_file_name, calibration_file_path),
                         "The data set should be reloaded")

        # Clean up
        delete_data_set(calibration_file_name)

if __name__ == "__main__":
    unittest.main()
