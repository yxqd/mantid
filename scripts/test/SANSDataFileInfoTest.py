import unittest

import mantid
from SANSDataFileInfo import SANSScatteringFileInfo


SANS_DATA_FILE_INFO_SAMPLE_DATA_SET_NAME = "test_workspace_sans_data_file_info"


class TestSANSScatteringFileInfo(unittest.TestCase):
    def _setup_default_info(self, has_dark_run=False, is_valid=True):
        info = SANSScatteringFileInfo()

        info.workspace_name = "test_name"
        info.history = "test_history"
        info.size = 123

        info.has_dark_run = has_dark_run

        info.tube_calibration_workspace_name = "test_tube_name"
        info.tube_calibration_file_name = "test_tube_file_name"

        info.geometry_shape = "test_shape"
        info.geometry_height = 123.3
        info.geometry_thickness = 23.4
        info.geometry_width = 1234.5

        info.beam_centre_position_1 = 123.4
        if is_valid:
            info.beam_centre_position_2 = 34.5

        return info

    def test_that_properties_raise_if_incorrect_type_is_passed(self):
        # Arrange
        info = SANSScatteringFileInfo()

        # Act + Assert
        self.assertRaises(ValueError, setattr, info, "workspace_name", 123)
        self.assertRaises(ValueError, setattr, info, "history", 123)
        self.assertRaises(ValueError, setattr, info, "size", "123")

        self.assertRaises(ValueError, setattr, info, "has_dark_run", "123")

        self.assertRaises(ValueError, setattr, info, "tube_calibration_workspace_name", 123)
        self.assertRaises(ValueError, setattr, info, "tube_calibration_file_name", 123)

        self.assertRaises(ValueError, setattr, info, "geometry_shape", 123)
        self.assertRaises(ValueError, setattr, info, "geometry_thickness", "123")
        self.assertRaises(ValueError, setattr, info, "geometry_width", "123")
        self.assertRaises(ValueError, setattr, info, "geometry_height", "123")

        self.assertRaises(ValueError, setattr, info, "beam_centre_position_1", "123")
        self.assertRaises(ValueError, setattr, info, "beam_centre_position_2", "123")

    def test_old_file_with_dark_run_enabled_requires_reload(self):
        # Arrange
        info_old = self._setup_default_info(has_dark_run=True, is_valid=True)
        info_new = self._setup_default_info(has_dark_run=False, is_valid=True)
        # Act
        requires_reload = info_old.requires_a_reload(info_new)
        # Assert
        self.assertTrue(requires_reload, "Having run a dark run should require a reload")

    def test_only_new_file_with_dark_run_enabled_does_no_trequire_reload(self):
        pass

    def test_that_invalid_state_raises(self):
        pass

    def test_that_differing_parameter_hash_requires_reload(self):
        pass

    def test_that_same_parameter_hash_does_not_require_reload(self):
        pass

"""
class TestSANSScatterinFileInfo(unittest.TestCase):
    def test_that_incomplete_original_file_info_raises(self):
        # Arrange
        ws_name = load_sample_data_set()
        ws_info = SANSScatterinFileInfo()
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
"""
if __name__ == "__main__":
    unittest.main()
