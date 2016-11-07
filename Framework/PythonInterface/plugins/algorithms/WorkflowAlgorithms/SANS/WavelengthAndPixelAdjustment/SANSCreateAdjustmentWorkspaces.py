# pylint: disable=invalid-name

""" SANSCreateAdjustmentWorkspaces algorithm creates workspaces for pixel adjustment
    , wavelength adjustment and pixel-and-wavelength adjustment workspaces.
"""

from mantid.kernel import (Direction, Property, PropertyManagerProperty)
from mantid.api import (DataProcessorAlgorithm, MatrixWorkspaceProperty, AlgorithmFactory, PropertyMode, Progress)

from SANS2.Common.SANSConstants import SANSConstants
from SANS2.Common.SANSFunctions import create_unmanaged_algorithm
from SANS2.State.SANSStateBase import create_deserialized_sans_state_from_property_manager


class SANSCreateAdjustmentWorkspaces(DataProcessorAlgorithm):
    def category(self):
        return 'SANS\\Adjust'

    def summary(self):
        return 'Calculates wavelength adjustment, pixel adjustment workspaces and wavelength-and-pixel ' \
               'adjustment workspaces.'

    def PyInit(self):
        # ---------------
        # INPUT
        # ---------------
        # State
        self.declareProperty(PropertyManagerProperty('SANSState'),
                             doc='A property manager which fulfills the SANSState contract.')

        # Input workspaces
        self.declareProperty(MatrixWorkspaceProperty('TransmissionWorkspace', '',
                                                     optional=PropertyMode.Optional, direction=Direction.Input),
                             doc='The transmission workspace.')
        self.declareProperty(MatrixWorkspaceProperty('DirectWorkspace', '',
                                                     optional=PropertyMode.Optional, direction=Direction.Input),
                             doc='The direct workspace.')
        self.declareProperty(MatrixWorkspaceProperty('MonitorWorkspace', '',
                                                     optional=PropertyMode.Optional, direction=Direction.Input),
                             doc='The scatter monitor workspace. This workspace only contains monitors.')
        self.declareProperty(MatrixWorkspaceProperty('SampleData', '',
                                                     optional=PropertyMode.Optional, direction=Direction.Input),
                             doc='A workspace cropped to the detector to be reduced (the SAME as the input to Q1D). '
                                 'This used to verify the solid angle. The workspace is not modified, just inspected.')

        # The component
        self.declareProperty('Component', '', direction=Direction.Input, doc='Component which is being investigated.')

        # Slice factor for monitor
        self.declareProperty('SliceEventFactor', 1.0, direction=Direction.Input, doc='The slice factor for the monitor '
                                                                                     'normalization. This factor is the'
                                                                                     ' one obtained from event '
                                                                                     'slicing.')

        # ---------------
        # Output
        # ---------------
        self.declareProperty(MatrixWorkspaceProperty('OutputWorkspaceWavelengthAdjustment', '',
                                                     direction=Direction.Output),
                             doc='The workspace for wavelength-based adjustments.')

        self.declareProperty(MatrixWorkspaceProperty('OutputWorkspacePixelAdjustment', '',
                                                     direction=Direction.Output),
                             doc='The workspace for wavelength-based adjustments.')

        self.declareProperty(MatrixWorkspaceProperty('OutputWorkspaceWavelengthAndPixelAdjustment', '',
                                                     direction=Direction.Output),
                             doc='The workspace for, both, wavelength- and pixel-based adjustments.')

    def PyExec(self):
        # Read the state
        state_property_manager = self.getProperty("SANSState").value
        state = create_deserialized_sans_state_from_property_manager(state_property_manager)

        # --------------------------------------
        # Get the monitor normalization workspace
        # --------------------------------------
        monitor_normalization_workspace = self._get_monitor_normalization_workspace(state)

        # --------------------------------------
        # Get the calculated transmission
        # --------------------------------------
        calculated_transmission_workspace = self._get_calculated_transmission_workspace(state)

        # --------------------------------------
        # Get the wide angle correction workspace
        # --------------------------------------
        wave_length_and_pixel_adjustment_workspace = self._get_wide_angle_correction_workspace(
                                                                   calculated_transmission_workspace)  # noqa

        # --------------------------------------------
        # Get the full wavelength and pixel adjustment
        # --------------------------------------------
        wave_length_adjustment_workspace, \
        pixel_length_adjustment_workspace = self._get_wavelength_and_pixel_adjustment_workspaces(state,
                                                                            monitor_normalization_workspace,  # noqa
                                                                            calculated_transmission_workspace)  # noqa

        if wave_length_adjustment_workspace:
            self.setProperty("OutputWorkspaceWavelengthAdjustment", wave_length_adjustment_workspace)
        if pixel_length_adjustment_workspace:
            self.getProperty("OutputWorkspacePixelAdjustment", pixel_length_adjustment_workspace)
        if wave_length_and_pixel_adjustment_workspace:
            self.getProperty("OutputWorkspaceWavelengthAndPixelAdjustment", wave_length_and_pixel_adjustment_workspace)

    def _get_wavelength_and_pixel_adjustment_workspaces(self, state,
                                                        monitor_normalization_workspace,
                                                        calculated_transmission_workspace,):
        component = self.getProperty("Component").value
        wave_pixel_adjustment_name = "SANSCreateWavelengthAndPixelAdjustment"
        wave_pixel_adjustment_options = {"State": state,
                                         "Transmission": calculated_transmission_workspace,
                                         "MonitorNormalization": monitor_normalization_workspace,
                                         "OutputWorkspaceWavelengthAdjustment": "dummy1",
                                         "OutputWorkspacePixelAdjustment": "dummy2",
                                         "Component": component}

        wave_pixel_adjustment_alg = create_unmanaged_algorithm(wave_pixel_adjustment_name,
                                                               *wave_pixel_adjustment_options)
        wave_pixel_adjustment_alg.execute()
        wavelength_out = wave_pixel_adjustment_alg.getProperty("OutputWorkspaceWavelengthAdjustment").value
        pixel_out = wave_pixel_adjustment_alg.getProperty("OutputWorkspacePixelAdjustment").value
        return wavelength_out, pixel_out

    def _get_monitor_normalization_workspace(self, state):
        """
        Gets the monitor normalization workspace via the SANSNormalizeToMonitor algorithm

        :param state: a SANSState object.
        :return: the normalization workspace.
        """
        monitor_workspace = self.getProperty("MonitorWorkspace").value
        scale_factor = self.getProperty("SliceEventFactor").value

        normalize_name = "SANSNormalizeMonitor"
        normalize_option = {SANSConstants.input_workspace: monitor_workspace,
                            SANSConstants.output_workspace: SANSConstants.dummy,
                            'SANSState': state,
                            "ScaleFactor": scale_factor}
        normalize_alg = create_unmanaged_algorithm(normalize_name, **normalize_option)
        normalize_alg.execute()
        return self.getProperty(SANSConstants.output_workspace).value

    def _get_calculated_transmission_workspace(self, state):
        """
        Creates the fitted transmission workspace.

        @param state: a SANSState object.
        @return: the transmission fit.
        """
        transmission_workspace = self.getProperty("TransmissionWorkspace").value
        direct_workspace = self.getProperty("DirectWorkspace").value

        transmission_name = "SANSCalculateTransmission"
        transmission_options = {"TransmissionWorkspace": transmission_workspace,
                                "DirectWorkspace": direct_workspace,
                                "SANSState": state,
                                SANSConstants.output_workspace: SANSConstants.dummy}
        transmission_alg = create_unmanaged_algorithm(transmission_name, **transmission_options)
        transmission_alg.execute()
        return transmission_alg.getProperty(SANSConstants.output_workspace).value

    def _get_wide_angle_correction_workspace(self, calculated_transmission_workspace):
        sample_data = self.getProperty("SampleData").value
        wide_angle_name = "SANSWideAngleCorrection"
        wide_angle_options = {"SampleData": sample_data,
                              "TransmissionData": calculated_transmission_workspace,
                              SANSConstants.output_workspace: SANSConstants.dummy}
        wide_angle_alg = create_unmanaged_algorithm(wide_angle_name, **wide_angle_options)
        wide_angle_alg.execute()
        return wide_angle_alg.getProperty(SANSConstants.output_workspace).value

    def validateInputs(self):
        errors = dict()
        # Check that the input can be converted into the right state object
        state_property_manager = self.getProperty("SANSState").value
        try:
            state = create_deserialized_sans_state_from_property_manager(state_property_manager)
            state.property_manager = state_property_manager
            state.validate()
        except ValueError as err:
            errors.update({"SANSCreateAdjustmentWorkspaces": str(err)})
        return errors


# Register algorithm with Mantid
AlgorithmFactory.subscribe(SANSCreateAdjustmentWorkspaces)
