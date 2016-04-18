import numbers
from mantid.simpleapi import *



class SANSScatteringFileInfo(object):
    def __init__(self):
        super(SANSScatteringFileInfo, self).__init__()

        # Workspace parameters
        self._workspace_name = None
        self._history = None
        self._size = None

        # Dark Run parameters - this is the only case when we strictly require a reload
        self._has_dark_run = None

        # Tube Calibration paramters
        self._tube_calibration_workspace_name = None
        self._tube_calibration_file_path = None

        # Geometry which is stored in sample
        self._geometry_shape = None
        self._geometry_thickness = None
        self._geometry_width = None
        self._geometry_height = None

        # Beam Centre
        self._beam_position_1 = None
        self._beam_position_2 = None

    @property
    def workspace_name(self):
        return self._workspace_name

    @workspace_name.setter
    def workspace_name(self, value):
        if not isinstance(value, basestring):
            raise ValueError("The workspace name must be a string")
        self._workspace_name = value

    @property
    def history(self):
        return self._history

    @history.setter
    def history(self, value):
        if not isinstance(value, basestring):
            raise ValueError("The history must be a string")
        self._history = value

    @property
    def size(self):
        return self._size

    @size.setter
    def size(self, value):
        if not isinstance(value, numbers.Integral):
            raise ValueError("The size must be an integer.")
        self._size = value

    @property
    def has_dark_run(self):
        return self._has_dark_run

    @has_dark_run.setter
    def has_dark_run(self, value):
        if not isinstance(value, bool):
            raise ValueError("The dark-run flag must be a bool.")
        self._has_dark_run = value

    @property
    def tube_calibration_workspace_name(self):
        return self._tube_calibration_workspace_name

    @tube_calibration_workspace_name.setter
    def tube_calibration_workspace_name(self, value):
        if not isinstance(value, basestring):
            raise ValueError("The calibration workspace name must be a string")
        self._tube_calibration_workspace_name = value

    @property
    def tube_calibration_file_name(self):
        return self._tube_calibration_workspace_name

    @tube_calibration_file_name.setter
    def tube_calibration_file_name(self, value):
        if not isinstance(value, basestring):
            raise ValueError("The calibration file path must be a string")
        self._tube_calibration_file_name = value

    @property
    def geometry_shape(self):
        return self._geometry_shape

    @geometry_shape.setter
    def geometry_shape(self, value):
        if not isinstance(value, basestring):
            raise ValueError("The geometry shape name must be a string")
        self._geometry_shape = value 

    @property
    def geometry_thickness(self):
        return self._geometry_thickness

    @geometry_thickness.setter
    def geometry_thickness(self, value):
        if not isinstance(value, numbers.Real):
            raise ValueError("The geometry thickness must be a real.")
        self._geometry_thickness = value 

    @property
    def geometry_width(self):
        return self._geometry_width

    @geometry_width.setter
    def geometry_width(self, value):
        if not isinstance(value, numbers.Real):
            raise ValueError("The geometry width must be a real.")
        self._geometry_width = value 

    @property
    def geometry_height(self):
        return self._geometry_height

    @geometry_height.setter
    def geometry_height(self, value):
        if not isinstance(value, numbers.Real):
            raise ValueError("The geometry width must be a real.")
        self._geometry_height = value 

    @property
    def beam_centre_position_1(self):
        return self._beam_position_1

    @beam_centre_position_1.setter
    def beam_centre_position_1(self, value):
        if not isinstance(value, numbers.Real):
            raise ValueError("The beam centre position 1 must be a real.")
        self._beam_position_1 = value 

    @property
    def beam_centre_position_2(self):
        return self._beam_position_2

    @beam_centre_position_2.setter
    def beam_centre_position_2(self, value):
        if not isinstance(value, numbers.Real):
            raise ValueError("The beam centre position 2 must be a real.")
        self._beam_position_2 = value 

    def is_in_valid_state(self):
        if (self.workspace_name is None or
            self.history is None or
            self.size is None or 
            self.has_dark_run is None or
            self.geometry_shape is None or
            self.geometry_thickness is None or
            self.geometry_width is None or
            self.geometry_height is None or
            self.beam_position_1 is None or
            self.beam_position_2 is None):
            return False
        else:
            return True

    def provide_parameter_hash(self):
        # Tube Calibration paramters
        tube_calibration_workspace_name = "" if self.tube_calibration_workspace_name is None else self.tube_calibration_workspace_name
        tube_calibration_file_path = "" if self.tube_calibration_workspace_name is None else self.tube_calibration_workspace_name

        param_hash = hash(tuple([self.workspace_name ,self.history, self.size, self.has_dark_run, 
                                tube_calibration_workspace_name, tube_calibration_file_path,
                                self.geometry_shape, self.geometry_thickness,
                                self.geometry_width, self.geometry_height,
                                self.beam_position_1, self.beam_position_2]))
        return param_hash

    def requires_a_reload(self, other_sans_data_file_info):
        if not self.is_in_valid_state() or not other_sans_data_file_info.is_in_valid_state():
            raise ValueError("The SANSScatterFileInfo is not in a valid state.")

        # If a dark run has been performed then require a reload
        if self.has_dark_run():
            return True

        # Compare parameter hash
        param_hash = self.provide_parameter_hash()
        param_hash_other = other_sans_data_file_info.provide_parameter_hash()

        if param_hash == param_hash_other:
            return False
        else:
            return True

    def _compare_work_space_parametrs(self, other_sans_data_file_info):
        return (self.workspace_name == other_sans_data_file_info.workspace_name and
                self.size == other_sans_data_file_info.size and
                self.history == other_sans_data_file_info.history)


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

    def set_sample_transmission(self, workspace_name):
        pass

    def set_sample_direct(self, workspace_name):
        pass

