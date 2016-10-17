# pylint: disable=invalid-name

""" SANSCreateAdjustmentWorkspaces algorithm creates workspaces for pixel adjustment
    , wavelength adjustment and pixel-and-wavelength adjustment workspaces.
"""

from mantid.kernel import (Direction, IntBoundedValidator, FloatBoundedValidator, StringListValidator,
                           Property, PropertyManagerProperty)
from mantid.api import (DataProcessorAlgorithm, MatrixWorkspaceProperty, AlgorithmFactory, PropertyMode, Progress)

from SANS2.Common.SANSConstants import SANSConstants
from SANS2.Common.SANSEnumerations import (RebinType, RangeStepType)
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
        calculated_transmission_workspace = self._get_calculated_transmission_workspace(adjustment_state)

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
        monitor_workspace = self.getProperty("MonitorWorkspace").value
        scale_factor = self.getProperty("SliceEventFactor").value

        incident_monitor = adjustment_state.incident_monitor

        # Prompt peak settings
        use_prompt_peak_correction = adjustment_state.use_prompt_peak_correction
        prompt_peak_correction_min = adjustment_state.prompt_peak_correction_min
        prompt_peak_correction_max = adjustment_state.prompt_peak_correction_max

        # Background settings
        background_tof_general_start = adjustment_state.background_TOF_general_start
        background_tof_general_stop = adjustment_state.background_TOF_general_stop
        background_tof_monitor_start = adjustment_state.background_TOF_general_start
        background_tof_monitor_stop = adjustment_state.background_TOF_general_stop

        # Wavelength rebin settings
        rebin_type = adjustment_state.rebin_type
        wavelength_low = adjustment_state.wavelength_low
        wavelength_high = adjustment_state.wavelength_high
        wavelength_step = adjustment_state.wavelength_step
        wavelength_step_type = adjustment_state.wavelength_step_type

        normalize_name = "SANSNormalizeMonitor"
        normalize_option = {SANSConstants.input_workspace: monitor_workspace,
                            SANSConstants.output_workspace: SANSConstants.dummy,
                            "IncidentMonitorSpectrumNumber": incident_monitor,
                            "ScaleFactor": scale_factor,
                            "WavelengthLow": wavelength_low,
                            "WavelengthHigh": wavelength_high,
                            "WavelengthStep": wavelength_step}

        # Select the correct wavelength range step type
        if wavelength_step_type is RangeStepType.Log:
            normalize_option.update({"WavelengthStepType": "LOG"})
        else:
            normalize_option.update({"WavelengthStepType": "LIN"})

        # Select the correct rebin mode
        if rebin_type is RebinType.InterpolatingRebin:
            normalize_option.update({"RebinMode": "InterpolatingRebin"})
        else:
            normalize_option.update({"RebinMode": "Rebin"})

        # Add the prompt peak correction information
        if use_prompt_peak_correction:
            normalize_option.update({"PromptPeakCorrectionStart", prompt_peak_correction_min,
                                    "PromptPeakCorrectionStop", prompt_peak_correction_max})

        # Add background correction information
        if background_tof_monitor_start is not None \
                and background_tof_monitor_stop is not None\
                and incident_monitor in background_tof_monitor_start \
                and  incident_monitor in background_tof_monitor_stop:
            background_start = background_tof_monitor_start[incident_monitor]
            background_stop = background_tof_monitor_stop[incident_monitor]
        else:
            background_start = background_tof_general_start
            background_stop = background_tof_general_stop
        normalize_option.update({"FlatBackgroundCorrectionStart": background_start,
                                 "FlatBackgroundCorrectionStop": background_stop})

        normalize_alg = create_unmanaged_algorithm(normalize_name, **normalize_option)
        normalize_alg.execute()
        return self.getProperty(SANSConstants.output_workspace).value

    def _get_calculated_transmission_workspace(self, adjustment_state):
        transmission_workspace = self.getProperty("TransmissionWorkspace").value
        direct_workspace = self.getProperty("DirectWorkspace").value

        # Incident monitor
        incident_monitor_for_transmission_calculation = adjustment_state.incident_monitor_for_transmission_calculation

        # ROI settings - Note that either ROI or the transmission monitor is expected to be used. The ROI settings
        # take precedence over the transmission monitor
        transmission_radius_on_detector = adjustment_state.transmission_radius_on_detector
        transmission_roi_files = adjustment_state.transmission_roi_files
        transmission_mask_files = adjustment_state.transmission_mask_files

        # Transmission monitor settings
        default_transmission_monitor = adjustment_state.default_transmission_monitor
        transmission_monitor = adjustment_state.transmission_monitor

        # Prompt peak settings
        use_prompt_peak_correction = adjustment_state.use_prompt_peak_correction
        prompt_peak_correction_min = adjustment_state.prompt_peak_correction_min
        prompt_peak_correction_max = adjustment_state.prompt_peak_correction_max

        # Background settings
        background_tof_general_start = adjustment_state.background_TOF_general_start
        background_tof_general_stop = adjustment_state.background_TOF_general_stop
        background_tof_monitor_start = adjustment_state.background_TOF_general_start
        background_tof_monitor_stop = adjustment_state.background_TOF_general_stop
        background_tof_ROI_start = adjustment_state.background_TOF_ROI_start
        background_tof_ROI_stop = adjustment_state.background_TOF_ROI_stop

        # Wavelength rebin settings
        rebin_type = adjustment_state.rebin_type
        wavelength_low = adjustment_state.wavelength_low
        wavelength_high = adjustment_state.wavelength_high
        wavelength_step = adjustment_state.wavelength_step
        wavelength_step_type = adjustment_state.wavelength_step_type

        use_full_wavelength_range = adjustment_state.use_full_wavelength_range
        wavelength_range_low = adjustment_state.wavelength_range_low
        wavelength_range_high = adjustment_state.wavelength_range_high
        default_wavelength_range_low = adjustment_state.default_wavelength_range_low
        default_wavelength_range_high = adjustment_state.default_wavelength_range_high

        # Fit settings
        fit_type = adjustment_state.fit_type
        polynomial_order = adjustment_state.polynomial_order

        transmission_name = "SANSCalculateTransmission"
        transmission_options = {"TransmissionWorkspace": transmission_workspace,
                                "DirectWorkspace": direct_workspace,
                                SANSConstants.output_workspace: SANSConstants.dummy,
                                "IncidentMonitorSpectrumNumber": incident_monitor_for_transmission_calculation}

        if prompt_peak_correction_min and prompt_peak_correction_max:
            transmission_options.update({"PromptPeakCorrectionStart": prompt_peak_correction_min},
                                        {"PromptPeakCorrectionStop": prompt_peak_correction_max})

        if transmission_radius_on_detector is not None:
            transmission_options.update({"TransmissionRadius": transmission_radius_on_detector})
        if transmission_roi_files:
            transmission_options.update({"TransmissionROIFiles": transmission_roi_files})
        if transmission_mask_files:
            transmission_options.update({"TransmissionMaskFiles": transmission_mask_files})

        if transmission_monitor is None:
            transmission_monitor = default_transmission_monitor
        transmission_options.update({"TransmissionMaskFiles": transmission_monitor})

        if background_tof_general_start and background_tof_general_stop:
            transmission_options.update({""})



    def validateInputs(self):
        errors = dict()

        return errors

# Register algorithm with Mantid
AlgorithmFactory.subscribe(SANSCreateAdjustmentWorkspaces)
