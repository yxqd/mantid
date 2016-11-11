import unittest
import mantid

from mantid.kernel import (PropertyManagerProperty, PropertyManager)
from mantid.api import Algorithm

from SANS2.State.SANSStateAdjustment import (SANSStateAdjustment, SANSStateAdjustmentISIS)
from SANS2.State.SANSStateBase import create_deserialized_sans_state_from_property_manager
from SANS2.State.StateDirector.TestDirector import TestDirector


class SANSStateReductionTest(unittest.TestCase):
    def test_that_is_sans_state_data_object(self):
        state = SANSStateAdjustmentISIS()
        self.assertTrue(isinstance(state, SANSStateAdjustment))

    def test_that_reduction_state_gets_and_sets(self):
        # Arrange
        state = SANSStateAdjustmentISIS()

        # Act
        # Get sample states from the test director, else we have to setup the sub states by hand.
        test_director = TestDirector()
        sample_state = test_director.construct()

        state.wide_angle_correction = True
        state.calculate_transmission = sample_state.adjustment.calculate_transmission
        state.normalize_to_monitor = sample_state.adjustment.normalize_to_monitor
        state.wavelength_and_pixel_adjustment = sample_state.adjustment.wavelength_and_pixel_adjustment

        # Assert
        self.assertTrue(state.wide_angle_correction)
        # The sub states are not tested since it is enough to call the validate method below to verify that they have
        # been set.

        try:
            state.validate()
            is_valid = True
        except ValueError:
            is_valid = False
        self.assertTrue(is_valid)

    def test_that_invalid_types_for_parameters_raise_type_error(self):
        # Arrange
        state = SANSStateAdjustmentISIS()

        # Act and Assert
        try:
            state.wide_angle_correction = ["sdf"]
            is_valid = True
        except TypeError:
            is_valid = False
        self.assertFalse(is_valid)

    def test_that_dict_can_be_generated_from_state_object_and_property_manager_read_in(self):
        class FakeAlgorithm(Algorithm):
            def PyInit(self):
                self.declareProperty(PropertyManagerProperty("Args"))

            def PyExec(self):
                pass

        # Arrange
        state = SANSStateAdjustmentISIS()

        # Get sample states from the test director, else we have to setup the sub states by hand.
        test_director = TestDirector()
        sample_state = test_director.construct()

        state.wide_angle_correction = False
        state.calculate_transmission = sample_state.adjustment.calculate_transmission
        state.normalize_to_monitor = sample_state.adjustment.normalize_to_monitor
        state.wavelength_and_pixel_adjustment = sample_state.adjustment.wavelength_and_pixel_adjustment

        # Act
        serialized = state.property_manager
        fake = FakeAlgorithm()
        fake.initialize()
        fake.setProperty("Args", serialized)
        property_manager = fake.getProperty("Args").value

        # Assert
        self.assertTrue(type(serialized) == dict)
        self.assertTrue(type(property_manager) == PropertyManager)
        state_2 = create_deserialized_sans_state_from_property_manager(property_manager)
        state_2.property_manager = property_manager

        self.assertTrue(not state_2.wide_angle_correction)
        try:
            state_2.validate()
            is_valid = True
        except ValueError:
            is_valid = False
        self.assertTrue(is_valid)


if __name__ == '__main__':
    unittest.main()
