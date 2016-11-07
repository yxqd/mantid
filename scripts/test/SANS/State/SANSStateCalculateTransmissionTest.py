import unittest
import mantid
from mantid.kernel import (PropertyManagerProperty, PropertyManager)
from mantid.api import Algorithm
from SANS2.State.SANSStateCalculateTransmission import (SANSStateCalculateTransmission,
                                                        SANSStateCalculateTransmissionISIS,
                                                        SANSStateCalculateTransmissionLOQ)
from SANS2.Common.SANSEnumerations import (RebinType, RangeStepType, FitType, DataType,
                                           convert_reduction_data_type_to_string)


class SANSStateCalculateTransmissionTest(unittest.TestCase):
    def test_that_is_sans_state_data_object(self):
        state = SANSStateCalculateTransmissionISIS()
        self.assertTrue(isinstance(state, SANSStateCalculateTransmission))

    def test_that_can_set_and_get_values(self):
        # Arrange
        state = SANSStateCalculateTransmissionISIS()

        # Act + Assert
        state.prompt_peak_correction_min = 12.0
        state.prompt_peak_correction_max = 17.0

        state.rebin_type = RebinType.Rebin
        state.wavelength_low = 1.5
        state.wavelength_high = 2.7
        state.wavelength_step = 0.5
        state.wavelength_step_type = RangeStepType.Lin
        state.use_full_wavelength_range = True
        state.wavelength_full_range_low = 123.0
        state.wavelength_full_range_high = 343.0

        state.incident_monitor = 1
        state.transmission_monitor = 2
        state.default_transmission_monitor = 3
        state.roi_files = ["sfsdf", "sdfsdfsf"]
        state.transmission_radius = 2.6

        state.background_TOF_general_start = 1.4
        state.background_TOF_general_stop = 34.6
        state.background_TOF_monitor_start = {"1": 123, "2": 123}
        state.background_TOF_monitor_stop = {"1": 234, "2": 2323}
        state.background_TOF_roi_start = 3.4
        state.background_TOF_roi_stop = 6.8

        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].wavelength_low = 10.0
        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].wavelength_high = 20.0
        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].polynomial_order = 3
        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].fit_type = FitType.Polynomial

        state.fit[convert_reduction_data_type_to_string(DataType.Can)].wavelength_low = 10.0
        state.fit[convert_reduction_data_type_to_string(DataType.Can)].wavelength_high = 20.0
        state.fit[convert_reduction_data_type_to_string(DataType.Can)].polynomial_order = 0
        state.fit[convert_reduction_data_type_to_string(DataType.Can)].fit_type = FitType.Linear

        try:
            state.validate()
            is_valid = True
        except ValueError:
            is_valid = False
        self.assertTrue(is_valid)

    def test_that_invalid_types_for_parameters_raise_type_error(self):
        # Arrange
        state = SANSStateCalculateTransmissionISIS()

        # Act + Assert
        try:
            state.background_TOF_monitor_start = "w234234"
            is_valid = True
        except TypeError:
            is_valid = False
        self.assertFalse(is_valid)

    def test_that_invalid_list_values_raise_value_error(self):
        # Arrange
        state = SANSStateCalculateTransmissionISIS()

        # Act + Assert
        try:
            state.wavelength_low = -1.0
            is_valid = True
        except ValueError:
            is_valid = False
        self.assertFalse(is_valid)

    def test_validate_method_raises_value_error_for_mismatching_monitor_start_and_stop_backgrounds(self):
        # Arrange
        state = SANSStateCalculateTransmissionISIS()

        state.prompt_peak_correction_min = 12.0
        state.prompt_peak_correction_max = 17.0

        state.rebin_type = RebinType.Rebin
        state.wavelength_low = 1.5
        state.wavelength_high = 2.7
        state.wavelength_step = 0.5
        state.wavelength_step_type = RangeStepType.Lin
        state.use_full_wavelength_range = True
        state.wavelength_full_range_low = 123.0
        state.wavelength_full_range_high = 343.0

        state.incident_monitor = 1
        state.transmission_monitor = 2
        state.default_transmission_monitor = 3
        state.roi_files = ["sfsdf", "sdfsdfsf"]
        state.transmission_radius = 2.6

        state.background_TOF_general_start = 1.4
        state.background_TOF_general_stop = 34.6
        state.background_TOF_monitor_start = {"1": 12, "2": 13}
        state.background_TOF_monitor_stop = {"1": 14, "3": 16}
        state.background_TOF_roi_start = 3.4
        state.background_TOF_roi_stop = 6.8

        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].wavelength_low = 10.0
        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].wavelength_high = 20.0
        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].polynomial_order = 3
        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].fit_type = FitType.Polynomial

        state.fit[convert_reduction_data_type_to_string(DataType.Can)].wavelength_low = 10.0
        state.fit[convert_reduction_data_type_to_string(DataType.Can)].wavelength_high = 20.0
        state.fit[convert_reduction_data_type_to_string(DataType.Can)].polynomial_order = 0
        state.fit[convert_reduction_data_type_to_string(DataType.Can)].fit_type = FitType.Linear

        # Act + Assert
        self.assertRaises(ValueError, state.validate)

    def test_validate_method_raises_when_no_transmission_is_specified(self):
        # Arrange
        state = SANSStateCalculateTransmissionISIS()

        state.prompt_peak_correction_min = 12.0
        state.prompt_peak_correction_max = 17.0

        state.rebin_type = RebinType.Rebin
        state.wavelength_low = 1.5
        state.wavelength_high = 2.7
        state.wavelength_step = 0.5
        state.wavelength_step_type = RangeStepType.Lin
        state.use_full_wavelength_range = True
        state.wavelength_full_range_low = 123.0
        state.wavelength_full_range_high = 343.0

        state.incident_monitor = 1

        state.background_TOF_general_start = 1.4
        state.background_TOF_general_stop = 34.6
        state.background_TOF_monitor_start = {"1": 12, "2": 13}
        state.background_TOF_monitor_stop = {"1": 14, "2": 16}
        state.background_TOF_roi_start = 3.4
        state.background_TOF_roi_stop = 6.8

        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].wavelength_low = 10.0
        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].wavelength_high = 20.0
        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].polynomial_order = 3
        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].fit_type = FitType.Polynomial

        state.fit[convert_reduction_data_type_to_string(DataType.Can)].wavelength_low = 10.0
        state.fit[convert_reduction_data_type_to_string(DataType.Can)].wavelength_high = 20.0
        state.fit[convert_reduction_data_type_to_string(DataType.Can)].polynomial_order = 0
        state.fit[convert_reduction_data_type_to_string(DataType.Can)].fit_type = FitType.Linear

        # Act + Assert
        self.assertRaises(ValueError, state.validate)

    def test_validate_method_raises_value_error_if_start_is_larger_than_stop(self):
        # Arrange
        state = SANSStateCalculateTransmissionISIS()
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
        state = SANSStateCalculateTransmissionLOQ()

        state.rebin_type = RebinType.Rebin
        state.wavelength_low = 1.5
        state.wavelength_high = 2.7
        state.wavelength_step = 0.5
        state.wavelength_step_type = RangeStepType.Lin
        state.use_full_wavelength_range = True
        state.wavelength_full_range_low = 123.0
        state.wavelength_full_range_high = 343.0

        state.incident_monitor = 1
        state.transmission_monitor = 2
        state.default_transmission_monitor = 3
        state.transmission_roi_files = ["sfsdf", "sdfsdfsf"]
        state.transmission_radius_on_detector = 2.6

        state.background_TOF_general_start = 1.4
        state.background_TOF_general_stop = 34.6
        state.background_TOF_monitor_start = {"1": 123, "2": 123}
        state.background_TOF_monitor_stop = {"1": 234, "2": 2323}
        state.background_TOF_roi_start = 3.4
        state.background_TOF_roi_stop = 6.8

        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].wavelength_low = 10.0
        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].wavelength_high = 20.0
        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].polynomial_order = 3
        state.fit[convert_reduction_data_type_to_string(DataType.Sample)].fit_type = FitType.Polynomial

        state.fit[convert_reduction_data_type_to_string(DataType.Can)].wavelength_low = 10.0
        state.fit[convert_reduction_data_type_to_string(DataType.Can)].wavelength_high = 20.0
        state.fit[convert_reduction_data_type_to_string(DataType.Can)].polynomial_order = 0
        state.fit[convert_reduction_data_type_to_string(DataType.Can)].fit_type = FitType.Linear

        # Act
        serialized = state.property_manager
        fake = FakeAlgorithm()
        fake.initialize()
        fake.setProperty("Args", serialized)
        property_manager = fake.getProperty("Args").value

        # Assert
        self.assertTrue(type(serialized) == dict)
        self.assertTrue(type(property_manager) == PropertyManager)
        state_2 = SANSStateCalculateTransmissionISIS()
        state_2.property_manager = property_manager

        self.assertTrue(state_2.prompt_peak_correction_min == 19000.0)
        self.assertTrue(state_2.prompt_peak_correction_max == 20500.0)

        self.assertTrue(state_2.rebin_type is RebinType.Rebin)
        self.assertTrue(state_2.wavelength_low == 1.5)
        self.assertTrue(state_2.wavelength_high == 2.7)
        self.assertTrue(state_2.wavelength_step == 0.5)
        self.assertTrue(state_2.wavelength_step_type is RangeStepType.Lin)
        self.assertTrue(state_2.use_full_wavelength_range is True)
        self.assertTrue(state_2.wavelength_full_range_low == 123.0)
        self.assertTrue(state_2.wavelength_full_range_high == 343.0)

        self.assertTrue(state_2.incident_monitor == 1)
        self.assertTrue(state_2.transmission_monitor == 2)
        self.assertTrue(state_2.default_transmission_monitor == 3)
        self.assertTrue(state_2.transmission_roi_files == ["sfsdf", "sdfsdfsf"])
        self.assertTrue(state_2.transmission_radius_on_detector == 2.6)

        self.assertTrue(state_2.background_TOF_general_start == 1.4)
        self.assertTrue(state_2.background_TOF_general_stop == 34.6)
        self.assertTrue(state_2.background_TOF_roi_start == 3.4)
        self.assertTrue(state_2.background_TOF_roi_stop == 6.8)

        self.assertTrue(len(set(state_2.background_TOF_monitor_start.items()) & set({"1": 123, "2": 123}.items())) == 2)
        self.assertTrue(len(set(state_2.background_TOF_monitor_stop.items()) & set({"1": 234, "2": 2323}.items())) == 2)

        self.assertTrue(state.fit[convert_reduction_data_type_to_string(DataType.Sample)].wavelength_low == 10.0)
        self.assertTrue(state.fit[convert_reduction_data_type_to_string(DataType.Sample)].wavelength_high == 20.0)
        self.assertTrue(state.fit[convert_reduction_data_type_to_string(DataType.Sample)].polynomial_order == 3)
        self.assertTrue(state.fit[convert_reduction_data_type_to_string(DataType.Sample)].fit_type
                        is FitType.Polynomial)

        self.assertTrue(state.fit[convert_reduction_data_type_to_string(DataType.Can)].wavelength_low == 10.0)
        self.assertTrue(state.fit[convert_reduction_data_type_to_string(DataType.Can)].wavelength_high == 20.0)
        self.assertTrue(state.fit[convert_reduction_data_type_to_string(DataType.Can)].polynomial_order == 0)
        self.assertTrue(state.fit[convert_reduction_data_type_to_string(DataType.Can)].fit_type is FitType.Linear)


if __name__ == '__main__':
    unittest.main()
