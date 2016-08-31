import unittest
import mantid
from mantid.kernel import (PropertyManagerProperty, PropertyManager)
from mantid.api import Algorithm
from SANS2.State.SANSStateScale import (SANSStateScale, SANSStateScaleISIS)
from SANS2.Common.SANSEnumerations import SampleShape


class SANSStateScaleTest(unittest.TestCase):
    def test_that_is_sans_state_data_object(self):
        state = SANSStateScaleISIS()
        self.assertTrue(isinstance(state, SANSStateScale))

    def test_that_can_set_and_get_values(self):
        # Arrange
        state = SANSStateScaleISIS()

        # Act + Assert
        state.shape = SampleShape.Cuboid
        state.thickness = 1.0
        state.width = 2.0
        state.height = 3.0
        state.scale = 4.0

        try:
            state.validate()
            is_valid = True
        except ValueError:
            is_valid = False
        self.assertTrue(is_valid)

    def test_that_invalid_types_for_parameters_raise_type_error(self):
        # Arrange
        state = SANSStateScaleISIS()

        # Act + Assert
        try:
            state.shape = "w234234"
            is_valid = True
        except TypeError:
            is_valid = False
        self.assertFalse(is_valid)

    def test_that_invalid_list_values_raise_value_error(self):
        # Arrange
        state = SANSStateScaleISIS()

        # Act + Assert
        try:
            state.thickness = -123.5
            is_valid = True
        except ValueError:
            is_valid = False
        self.assertFalse(is_valid)

    def test_that_dict_can_be_generated_from_state_object_and_property_manager_read_in(self):
        class FakeAlgorithm(Algorithm):
            def PyInit(self):
                self.declareProperty(PropertyManagerProperty("Args"))

            def PyExec(self):
                pass

        # Arrange
        state = SANSStateScaleISIS()
        shape = SampleShape.Cuboid
        thickness = 1.0
        height = 2.0
        width = 3.0
        scale = 4.0
        state.shape = shape
        state.scale = scale
        state.thickness = thickness
        state.height = height
        state.width = width

        # Act
        serialized = state.property_manager
        fake = FakeAlgorithm()
        fake.initialize()
        fake.setProperty("Args", serialized)
        property_manager = fake.getProperty("Args").value

        # Assert
        self.assertTrue(type(serialized) == dict)
        self.assertTrue(type(property_manager) == PropertyManager)
        state_2 = SANSStateScaleISIS()
        state_2.property_manager = property_manager

        self.assertTrue(state_2.shape is SampleShape.Cuboid)
        self.assertTrue(state_2.scale == scale)
        self.assertTrue(state_2.height == height)
        self.assertTrue(state_2.thickness == thickness)
        self.assertTrue(state_2.width == width)

if __name__ == '__main__':
    unittest.main()
