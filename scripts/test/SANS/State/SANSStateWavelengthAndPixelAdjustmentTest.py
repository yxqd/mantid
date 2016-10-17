import unittest
import mantid
from mantid.kernel import (PropertyManagerProperty, PropertyManager)
from mantid.api import Algorithm
from SANS2.State.SANSStateWavelengthAndPixelAdjustment import (SANSStateWavelengthAndPixelAdjustment,
                                                               SANSStateWavelengthAndPixelAdjustmentISIS)
from SANS2.Common.SANSEnumerations import (RebinType, RangeStepType)


class SANSStateWavelengthAndPixelAdjustmentTest(unittest.TestCase):
    def test_that_is_sans_state_data_object(self):
        state = SANSStateWavelengthAndPixelAdjustmentISIS()
        self.assertTrue(isinstance(state, SANSStateWavelengthAndPixelAdjustment))

    def test_that_can_set_and_get_values(self):
        # Arrange
        state = SANSStateWavelengthAndPixelAdjustmentISIS()

        # Act + Assert
        state.pixel_adjustment_file = "tests"
        state.wavelength_adjustment_file = "tests2"

        state.wavelength_low = 1.5
        state.wavelength_high = 2.7
        state.wavelength_step = 0.5
        state.wavelength_step_type = RangeStepType.Lin

        try:
            state.validate()
            is_valid = True
        except ValueError:
            is_valid = False
        self.assertTrue(is_valid)

    def test_that_invalid_types_for_parameters_raise_type_error(self):
        # Arrange
        state = SANSStateWavelengthAndPixelAdjustmentISIS()

        # Act + Assert
        try:
            state.pixel_adjustment_file = 12.5
            is_valid = True
        except TypeError:
            is_valid = False
        self.assertFalse(is_valid)

    def test_that_invalid_list_values_raise_value_error(self):
        # Arrange
        state = SANSStateWavelengthAndPixelAdjustmentISIS()

        # Act + Assert
        try:
            state.wavelength_low = -1.0
            is_valid = True
        except ValueError:
            is_valid = False
        self.assertFalse(is_valid)

    def test_validate_method_raises_value_error_if_start_is_larger_than_stop(self):
        # Arrange
        state = SANSStateWavelengthAndPixelAdjustmentISIS()
        state.wavelength_low = 1.5
        state.wavelength_high = 0.1
        state.wavelength_step = 0.5
        state.wavelength_step_type = RangeStepType.Lin

        # Act + Assert
        self.assertRaises(ValueError, state.validate)

    def test_that_dict_can_be_generated_from_state_object_and_property_manager_read_in(self):
        class FakeAlgorithm(Algorithm):
            def PyInit(self):
                self.declareProperty(PropertyManagerProperty("Args"))

            def PyExec(self):
                pass

        # Arrange
        state = SANSStateWavelengthAndPixelAdjustmentISIS()
        state.pixel_adjustment_file = "tests"
        state.wavelength_adjustment_file = "tests2"

        state.wavelength_low = 1.5
        state.wavelength_high = 2.7
        state.wavelength_step = 0.5
        state.wavelength_step_type = RangeStepType.Lin

        # Act
        serialized = state.property_manager
        fake = FakeAlgorithm()
        fake.initialize()
        fake.setProperty("Args", serialized)
        property_manager = fake.getProperty("Args").value

        # Assert
        self.assertTrue(type(serialized) == dict)
        self.assertTrue(type(property_manager) == PropertyManager)
        state_2 = SANSStateWavelengthAndPixelAdjustmentISIS()
        state_2.property_manager = property_manager

        self.assertTrue(state_2.pixel_adjustment_file == "tests")
        self.assertTrue(state_2.wavelength_adjustment_file == "tests2")
        self.assertTrue(state_2.wavelength_step_type == RangeStepType.Lin)
        self.assertTrue(state_2.wavelength_low == 1.5)
        self.assertTrue(state_2.wavelength_high == 2.7)
        self.assertTrue(state_2.wavelength_step == 0.5)
        self.assertTrue(state_2.wavelength_low == 1.5)


if __name__ == '__main__':
    unittest.main()
