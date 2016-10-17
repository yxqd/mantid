import unittest
import mantid

from SANS2.State.StateBuilder.SANSStateDataBuilder import get_data_builder
from SANS2.State.StateBuilder.SANSStateWavelengthAndPixelAdjustmentBuilder import get_wavelength_and_pixel_adjustment_builder
from SANS2.Common.SANSEnumerations import (RebinType, RangeStepType)
from SANS2.Common.SANSEnumerations import (SANSFacility, SANSInstrument)


class SANSStateReductionBuilderTest(unittest.TestCase):
    def test_that_reduction_state_can_be_built(self):
        # Arrange
        facility = SANSFacility.ISIS
        data_builder = get_data_builder(facility)
        data_builder.set_sample_scatter("LOQ74044")
        data_info = data_builder.build()

        # Act
        builder = get_wavelength_and_pixel_adjustment_builder(data_info)
        self.assertTrue(builder)

        builder.set_pixel_adjustment_file("test")
        builder.set_wavelength_adjustment_file("testst")
        builder.set_wavelength_low(1.5)
        builder.set_wavelength_high(2.7)
        builder.set_wavelength_step(0.5)
        builder.set_wavelength_step_type(RangeStepType.Lin)

        state = builder.build()

        # Assert
        self.assertTrue(state.pixel_adjustment_file == "test")
        self.assertTrue(state.wavelength_adjustment_file == "testst")
        self.assertTrue(state.wavelength_low == 1.5)
        self.assertTrue(state.wavelength_high == 2.7)
        self.assertTrue(state.wavelength_step == 0.5)
        self.assertTrue(state.wavelength_step_type is RangeStepType.Lin)


if __name__ == '__main__':
    unittest.main()
