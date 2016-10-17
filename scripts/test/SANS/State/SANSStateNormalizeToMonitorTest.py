import unittest
import mantid
from mantid.kernel import (PropertyManagerProperty, PropertyManager)
from mantid.api import Algorithm
from SANS2.State.SANSStateNormalizeToMonitor import (SANSStateNormalizeToMonitor, SANSStateNormalizeToMonitorISIS,
                                                     SANSStateNormalizeToMonitorLOQ)
from SANS2.Common.SANSEnumerations import (RebinType, RangeStepType)


class SANSStateNormalizeToMonitorTest(unittest.TestCase):
    def test_that_is_sans_state_data_object(self):
        state = SANSStateNormalizeToMonitorISIS()
        self.assertTrue(isinstance(state, SANSStateNormalizeToMonitor))

    def test_that_can_set_and_get_values(self):
        # Arrange
        state = SANSStateNormalizeToMonitorISIS()

        # Act + Assert
        state.prompt_peak_correction_min = 12.0
        state.prompt_peak_correction_max = 17.0

        state.rebin_type = RebinType.Rebin
        state.wavelength_low = 1.5
        state.wavelength_high = 2.7
        state.wavelength_step = 0.5
        state.wavelength_step_type = RangeStepType.Lin
        state.incident_monitor = 1

        state.background_TOF_general_start = 1.4
        state.background_TOF_general_stop = 34.6
        state.background_TOF_monitor_start = {"1": 123, "2": 123}
        state.background_TOF_monitor_stop = {"1": 234, "2": 2323}

        try:
            state.validate()
            is_valid = True
        except ValueError:
            is_valid = False
        self.assertTrue(is_valid)

    def test_that_invalid_types_for_parameters_raise_type_error(self):
        # Arrange
        state = SANSStateNormalizeToMonitorISIS()

        # Act + Assert
        try:
            state.background_TOF_monitor_start = "w234234"
            is_valid = True
        except TypeError:
            is_valid = False
        self.assertFalse(is_valid)

    def test_that_invalid_list_values_raise_value_error(self):
        # Arrange
        state = SANSStateNormalizeToMonitorISIS()

        # Act + Assert
        try:
            state.wavelength_low = -1.0
            is_valid = True
        except ValueError:
            is_valid = False
        self.assertFalse(is_valid)

    def test_validate_method_raises_value_error_for_mismatching_monitor_start_and_stop_backgrounds(self):
        # Arrange
        state = SANSStateNormalizeToMonitorISIS()
        state.background_TOF_monitor_start = {"1": 12, "2": 13}
        state.background_TOF_monitor_stop = {"1": 14, "3": 16}

        state.rebin_type = RebinType.Rebin
        state.wavelength_low = 1.5
        state.wavelength_high = 2.7
        state.wavelength_step = 0.5
        state.wavelength_step_type = RangeStepType.Lin
        state.incident_monitor = 1

        # Act + Assert
        self.assertRaises(ValueError, state.validate)

    def test_validate_method_raises_value_error_if_start_is_larger_than_stop(self):
        # Arrange
        state = SANSStateNormalizeToMonitorISIS()
        state.background_TOF_monitor_start = {"1": 12, "2": 13}
        state.background_TOF_monitor_stop = {"1": 11, "2": 13}
        state.rebin_type = RebinType.Rebin
        state.wavelength_low = 1.5
        state.wavelength_high = 2.7
        state.wavelength_step = 0.5
        state.wavelength_step_type = RangeStepType.Lin
        state.incident_monitor = 1

        # Act + Assert
        self.assertRaises(ValueError, state.validate)

    def test_that_dict_can_be_generated_from_state_object_and_property_manager_read_in(self):
        class FakeAlgorithm(Algorithm):
            def PyInit(self):
                self.declareProperty(PropertyManagerProperty("Args"))

            def PyExec(self):
                pass

        # Arrange
        state = SANSStateNormalizeToMonitorLOQ()

        state.background_TOF_monitor_start = {"1": 12, "2": 13}
        state.background_TOF_monitor_stop = {"1": 14, "2": 13}

        state.rebin_type = RebinType.Rebin
        state.wavelength_low = 1.5
        state.wavelength_high = 2.7
        state.wavelength_step = 0.5
        state.wavelength_step_type = RangeStepType.Lin
        state.incident_monitor = 1

        # Act
        serialized = state.property_manager
        fake = FakeAlgorithm()
        fake.initialize()
        fake.setProperty("Args", serialized)
        property_manager = fake.getProperty("Args").value

        # Assert
        self.assertTrue(type(serialized) == dict)
        self.assertTrue(type(property_manager) == PropertyManager)
        state_2 = SANSStateNormalizeToMonitorISIS()
        state_2.property_manager = property_manager

        self.assertTrue(state_2.rebin_type is RebinType.Rebin)
        self.assertTrue(state_2.wavelength_step_type == RangeStepType.Lin)
        self.assertTrue(state_2.wavelength_low == 1.5)
        self.assertTrue(state_2.wavelength_high == 2.7)
        self.assertTrue(state_2.wavelength_step == 0.5)
        self.assertTrue(state_2.wavelength_low == 1.5)
        self.assertTrue(state_2.incident_monitor == 1)

        # Checking the default value settings of LOQ
        self.assertTrue(state_2.prompt_peak_correction_min == 19000.0)
        self.assertTrue(state_2.prompt_peak_correction_max == 20500.0)

        start_2 = state_2.background_TOF_monitor_start
        stop_2 = state_2.background_TOF_monitor_stop

        self.assertTrue(len(set({"1": 12, "2": 13}.items()) & set(start_2.items())) == 2)
        self.assertTrue(len(set({"1": 14, "2": 13}.items()) & set(stop_2.items())) == 2)


if __name__ == '__main__':
    unittest.main()
