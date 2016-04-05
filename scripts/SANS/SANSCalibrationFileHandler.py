from mantid.simpleapi import *

class SANSCalibrationFileHandler(object):
    def __init__(self):
        super(SANSCalibrationFileHandler, self).__init__()
        self._calibration_file_name = None
        self._calibration_file_path = None
        self._calibration_was_updated= False

    @property
    def calibration_was_updated(self):
        return self._calibration_was_updated

    @calibration_was_updated.setter
    def calibration_was_updated(self, value):
        if isinstance(value, bool):
            self._calibration_was_updated = value
        else:
            raise ValueError("The calibration-was-updated flag requires a boolean value.")

    def provide_calibration_file(self, calibration_file_name, calibration_file_path):
        if self.requires_no_reload(calibration_file_name, calibration_file_path):
            # We can reuse the existing calibration file
            return mtd[self._calibration_file_name]
        else:
            # Need to load the calibration file fresh, reset the name and the upate behaviour
            self._calibration_file_name = calibration_file_name
            self._calibration_file_path = calibration_file_path
            self.calibration_was_updated = False
            return Load(Filename=calibration_file_path, OutputWorkspace=calibration_file_name)

    def requires_no_reload(self, calibration_file_name, calibration_file_path):
        return (self._calibration_file_exists_in_analysis_data_service(calibration_file_name) and
            self._calibration_file_names_are_identical(calibration_file_name) and
            self._calibration_file_paths_are_identical(calibration_file_path) and
              not self.calibration_was_updated)

    def _calibration_file_exists_in_analysis_data_service(self, calibration_file_name):
        return mtd.doesExist(calibration_file_name)

    def _calibration_file_names_are_identical(self, calibration_file_name):
        return calibration_file_name == self._calibration_file_name

    def _calibration_file_paths_are_identical(self, calibration_file_path):
        return calibration_file_path == self._calibration_file_path
