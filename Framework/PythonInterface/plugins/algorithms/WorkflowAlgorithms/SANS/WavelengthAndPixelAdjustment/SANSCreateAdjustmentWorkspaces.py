# pylint: disable=invalid-name

""" SANSCreateAdjustmentWorkspaces algorithm creates workspaces for pixel adjustment
    , wavelength adjustment and pixel-and-wavelength adjustment workspaces.
"""

from mantid.kernel import (Direction, IntBoundedValidator, FloatBoundedValidator, StringListValidator,
                           Property, PropertyManagerProperty)
from mantid.api import (DataProcessorAlgorithm, MatrixWorkspaceProperty, AlgorithmFactory, PropertyMode, Progress)

from SANS2.Common.SANSConstants import SANSConstants
from SANS2.Common.SANSFunctions import create_unmanaged_algorithm
from SANS2.State.SANSStateBase import create_deserialized_sans_state_from_property_manager


class SANSCreateAdjustmentWorkspaces(DataProcessorAlgorithm):
    def category(self):
        return 'SANS\\Normalize'

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

    def PyExec(self):
        # Read the state
        state_property_manager = self.getProperty("SANSState").value
        state = create_deserialized_sans_state_from_property_manager(state_property_manager)
        adjustment_state = state.adjustment
        # --------------------------------------
        # Get the monitor normalization workspace
        # --------------------------------------
        monitor_normalization_workspace = self._get_monitor_normalization_workspace(adjustment_state)

        # --------------------------------------
        # Get the calculated transmission
        # --------------------------------------
        calculated_transmission_workspace = self._get_calculated_transmission_workspace()

        # --------------------------------------
        # Get the wide angle correction workspace
        # --------------------------------------
        wide_angle_correction_workspace = self._get_wide_angle_correction_workspace()


    def _get_monitor_normalization_workspace(self, adjustment_state):
        """
        Gets the monitor normalization workspace via the SANSNormalizeToMonitor algorithm

        :param adjustment_state: a SANSStateAdjustment object.
        :return: the normalization workspace.
        """
        incident_monitor = adjustment_state.incident_monitor
        use_prompt_peak_correction = adjustment_state.use_prompt_peak_correction
        prompt_peak_correction_min = adjustment_state.prompt_peak_correction_min
        prompt_peak_correction_max = adjustment_state.prompt_peak_correction_max

        rebin_type = adjustment_state.rebin_type
        wavelength_low = adjustment_state.wavelength_low
        wavelength_high = adjustment_state.wavelength_high
        wavelength_step = adjustment_state.wavelength_step
        wavelength_step_type = adjustment_state.wavelength_step_type

        monitor_workspace = self.getProperty("MonitorWorkspace").value
        normalize_name = "SANSNormalizeMonitor"
        normalize_option = {SANSConstants.input_workspace}

    def validateInputs(self):
        errors = dict()

        return errors

# Register algorithm with Mantid
AlgorithmFactory.subscribe(SANSCreateWavelengthAndPixelAdjustment)
