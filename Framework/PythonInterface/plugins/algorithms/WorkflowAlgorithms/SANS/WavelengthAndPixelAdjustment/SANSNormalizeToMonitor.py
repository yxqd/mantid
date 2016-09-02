# pylint: disable=invalid-name

""" SANSNormalizeToMonitor algorithm calculates the normalization to the monitor."""

from mantid.kernel import (Direction, IntBoundedValidator, FloatBoundedValidator, StringListValidator, Property)
from mantid.api import (DataProcessorAlgorithm, MatrixWorkspaceProperty, AlgorithmFactory, PropertyMode,
                        FileProperty, FileAction, Progress)
from SANS2.Common.SANSConstants import SANSConstants
from SANS2.Common.SANSFunctions import create_unmanaged_algorithm


class SANSNormalizeToMonitor(DataProcessorAlgorithm):
    def category(self):
        return 'SANS\\Normalize'

    def summary(self):
        return 'Calculates a monitor normalization workspace for a SANS reduction.'

    def PyInit(self):
        # Input workspace in TOF
        self.declareProperty(MatrixWorkspaceProperty(SANSConstants.input_workspace, '',
                                                     optional=PropertyMode.Mandatory, direction=Direction.Input),
                             doc='The monitor workspace in time-of-flight units.')

        # Output workspace
        self.declareProperty(MatrixWorkspaceProperty(SANSConstants.output_workspace, '', direction=Direction.Output),
                             doc='A monitor normalization workspace in units of wavelength.')

        # The incident monitor
        self.declareProperty('IncidentMonitorSpectrumNumber', defaultValue=Property.EMPTY_INT,
                             direction=Direction.Input,
                             doc='The spectrum number of the input monitor')

        # A scale factor which could come from event workspace slicing. If the actual data workspace was sliced,
        # then one needs to scale the monitor measurement proportionally. This input is intended for this matter
        self.declareProperty('ScaleFactor', defaultValue=1.0, direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='Optional scale factor for the input workspace.')

        # Prompt Peak correction. This is a peak on the monitor data which comes from a burst of fast neutrons.
        # This peak is normally removed and an interpolation is performed between the start and end value of that
        # region.
        self.declareProperty('PromptPeakCorrectionStart', defaultValue=Property.EMPTY_DBL,
                             direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='The start time of the prompt peak.')
        self.declareProperty('PromptPeakCorrectionStop', defaultValue=Property.EMPTY_DBL,
                             direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='The stop time of the prompt peak.')

        # Flat background settings
        self.declareProperty('FlatBackgroundCorrectionStart', defaultValue=Property.EMPTY_DBL,
                             direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='The start time of the flat background correction.')
        self.declareProperty('FlatBackgroundCorrectionStop', defaultValue=Property.EMPTY_DBL,
                             direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='The stop time of the flat background correction.')

        # The output will be returned in wavelength units. These properties are for the wavelength rebinning.
        self.declareProperty('WavelengthLow', defaultValue=Property.EMPTY_DBL, direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='The low value of the wavelength binning.')
        self.declareProperty('WavelengthHigh', defaultValue=Property.EMPTY_DBL, direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='The high value of the wavelength binning.')
        self.declareProperty('WavelengthStep', defaultValue=Property.EMPTY_DBL, direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='The step size of the wavelength binning.')
        allowed_step_types = StringListValidator(["LOG", "LIN"])
        self.declareProperty('WavelengthStepType', "LIN", validator=allowed_step_types, direction=Direction.Input,
                             doc='The step type for rebinning.')
        allowed_rebin_methods = StringListValidator(["Rebin", "InterpolatingRebin"])
        self.declareProperty("RebinMode", "Rebin", validator=allowed_rebin_methods, direction=Direction.Input,
                             doc="The method which is to be applied to the rebinning.")

    def PyExec(self):
        # 1. Extract the spectrum of the incident monitor
        incident_monitor_spectrum_number = self.getProperty("IncidentMonitorSpectrumNumber").value
        workspace = self._extract_monitor(incident_monitor_spectrum_number)

        # 2. Multiply the workspace by the specified scaling factor.
        scale_factor = self.getProperty("ScaleFactor").value
        if scale_factor != 1.0:
            workspace = self._scale(workspace, scale_factor)

        # 3. Remove the prompt peak (if it exists)
        workspace = self._perform_prompt_peak_correction(workspace)

        # 4. Perform a flat background correction
        workspace = self._perform_flat_background_correction(workspace)

        # 5. Convert to wavelength with the specified bin settings.
        workspace = self._convert_to_wavelength(workspace)
        self.setProperty(SANSConstants.output_workspace, workspace)

    def _scale(self, workspace, factor):
        """
        The incident monitor is scaled by a factor.

        When we work with sliced event data, then we need to slice the monitor data set accordingly. The monitor
        needs to be scaled by the slice factor which one gets when operating SANSSliceEvent. If this was not performed,
        then we would be comparing the full monitor data with only parts of the detector data.
        :param workspace: the workspace to scale.
        :param factor: the scaling factor.
        :return: a scaled workspace.
        """
        scale_name = "Scale"
        scale_options = {SANSConstants.input_workspace: workspace,
                         SANSConstants.output_workspace: SANSConstants.dummy,
                         "Factor": factor,
                         "Operation": "Multiply"}
        scale_alg = create_unmanaged_algorithm(scale_name, **scale_options)
        scale_alg.execute()
        return scale_alg.getProperty(SANSConstants.output_workspace).value

    def _extract_monitor(self, spectrum_number):
        """
        The extracts a single spectrum from the input workspace.

        We are only interested in the incident monitor here.
        :param spectrum_number: the spectrum number of the incident beam monitor.
        :return: a workspace which only contains the incident beam spectrum.
        """
        workspace = self.getProperty(SANSConstants.input_workspace).value
        workspace_index = workspace.getIndexFromSpectrumNumber(spectrum_number)
        extract_name = "ExtractSingleSpectrum"
        extract_options = {SANSConstants.input_workspace: workspace,
                           SANSConstants.output_workspace: SANSConstants.dummy,
                           "WorkspaceIndex": workspace_index}
        extract_alg = create_unmanaged_algorithm(extract_name, **extract_options)
        extract_alg.execute()
        return extract_alg.getProperty(SANSConstants.output_workspace).value

    def _perform_prompt_peak_correction(self, workspace):
        """
        Performs a prompt peak correction.

        A prompt peak can occur when very fast neutrons shoot through the measurement. This can happen when working
        with two time regimes. Prompt peaks are prominent peaks which stand out from usual data. They occur frequently
        on LOQ, but are now also a possibility on other instruments. We deal with them, by removing the data and
        interpolating between the edge data points. If the user does not specify a start and stop time for the
        prompt peak, then this correction is not performed.
        :param workspace: the workspace which is to be corrected.
        :return: the corrected workspace.
        """
        prompt_peak_correction_start = self.getProperty("PromptPeakCorrectionStart").value
        prompt_peak_correction_stop = self.getProperty("PromptPeakCorrectionStop").value

        # We perform only a prompt peak correction if the start and stop values of the bins we want to remove,
        # were explicitly set. Some instruments require it, others don't.
        if prompt_peak_correction_start != Property.EMPTY_DBL and prompt_peak_correction_stop != Property.EMPTY_DBL:
            remove_name = "RemoveBins"
            remove_options = {SANSConstants.input_workspace: workspace,
                              SANSConstants.output_workspace: SANSConstants.dummy,
                              "XMin": prompt_peak_correction_start,
                              "XMax": prompt_peak_correction_stop,
                              "Interpolation": "Linear"}
            remove_alg = create_unmanaged_algorithm(remove_name, **remove_options)
            remove_alg.execute()
            workspace = remove_alg.getProperty(SANSConstants.output_workspace).value
        return workspace

    def _perform_flat_background_correction(self, workspace):
        """
        Removes an offset from the monitor data.

        A certain region of the data set is selected which corresponds to only background data. This data is averaged
        which results in a mean background value which is subtracted from the data.
        :param workspace: the workspace which is to be corrected.
        :return: the corrected workspace.
        """
        start_tof = self.getProperty("FlatBackgroundCorrectionStart").value
        stop_tof = self.getProperty("FlatBackgroundCorrectionStop").value
        flat_name = "CalculateFlatBackground"
        flat_options = {SANSConstants.input_workspace: workspace,
                        SANSConstants.output_workspace: SANSConstants.dummy,
                        "StartX": start_tof,
                        "EndX": stop_tof,
                        "Mode": "Mean"}
        flat_alg = create_unmanaged_algorithm(flat_name, **flat_options)
        flat_alg.execute()
        return flat_alg.getProperty(SANSConstants.output_workspace).value

    def _convert_to_wavelength(self, workspace):
        """
        Converts the workspace from time-of-flight units to wavelength units

        :param workspace: a time-of-flight workspace.
        :return: a wavelength workspace.
        """
        wavelength_low = self.getProperty("WavelengthLow").value
        wavelength_high = self.getProperty("WavelengthHigh").value
        wavelength_step = self.getProperty("WavelengthStep").value
        wavelength_step_type = self.getProperty("WavelengthStepType").value
        wavelength_rebin_mode = self.getProperty("RebinMode").value

        convert_name = "SANSConvertToWavelength"
        convert_options = {SANSConstants.input_workspace: workspace,
                           SANSConstants.output_workspace: SANSConstants.dummy,
                           "WavelengthLow": wavelength_low,
                           "WavelengthHigh": wavelength_high,
                           "WavelengthStep": wavelength_step,
                           "WavelengthStepType": wavelength_step_type,
                           "RebinMode": wavelength_rebin_mode}

        convert_alg = create_unmanaged_algorithm(convert_name, **convert_options)
        convert_alg.execute()
        return convert_alg.getProperty(SANSConstants.output_workspace).value

    def validateInputs(self):
        errors = dict()

        # ----------------------
        # Incident Monitor
        # ----------------------
        incident_monitor_spectrum_number= self.getProperty("IncidentMonitorSpectrumNumber").value

        # Ensure that the spectrum_number exists on the input workspace
        input_workspace = self.getProperty(SANSConstants.input_workspace).value
        try:
            input_workspace.getIndexFromSpectrumNumber(incident_monitor_spectrum_number)
            spectrum_number_exists_on_workspace = True
        except RuntimeError:
            spectrum_number_exists_on_workspace = False
        if not spectrum_number_exists_on_workspace:
            errors.update({"IncidentMonitorSpectrumNumber": "The incident monitor spectrum number {0} does not seem to "
                                                            "exist on the workspace."
                          .format(incident_monitor_spectrum_number)})

        # Ensure that the selected spectrum is actually a monitor
        if spectrum_number_exists_on_workspace:
            workspace_index = input_workspace.getIndexFromSpectrumNumber(incident_monitor_spectrum_number)
            detector = input_workspace.getDetector(workspace_index)
            if not detector.isMonitor():
                errors.update({"IncidentMonitorSpectrumNumber": "The incident monitor spectrum number does "
                                                                "not correspond to a monitor."})

        # ----------------------
        # Prompt Peak
        # ----------------------
        prompt_peak_correction_start = self.getProperty("PromptPeakCorrectionStart").value
        prompt_peak_correction_stop = self.getProperty("PromptPeakCorrectionStop").value

        # Either both prompt peak values are specified or they are not
        if (prompt_peak_correction_start == Property.EMPTY_DBL and prompt_peak_correction_stop != Property.EMPTY_DBL) \
                or (prompt_peak_correction_start != Property.EMPTY_DBL and
                    prompt_peak_correction_stop == Property.EMPTY_DBL):
            errors.update({"PromptPeakCorrectionStart": "Either both prompt peak correction entries "
                                                        "are specified or none."})
            errors.update({"PromptPeakCorrectionStop": "Either both prompt peak correction entries "
                                                       "are specified or none."})

        # The prompt peak start value needs to be smaller than the prompt peak stop value
        if (prompt_peak_correction_start != Property.EMPTY_DBL) and (prompt_peak_correction_stop != Property.EMPTY_DBL)\
                and prompt_peak_correction_start > prompt_peak_correction_stop:
            errors.update({"PromptPeakCorrectionStart": "The prompt peak start value needs to be smaller "
                                                        "than the stop value."})
            errors.update({"PromptPeakCorrectionStop": "The prompt peak start value needs to be smaller "
                                                       "than the stop value."})

        # --------------------------
        # Flat background correction
        # --------------------------
        flat_background_correction_start = self.getProperty("FlatBackgroundCorrectionStart").value
        flat_background_correction_stop = self.getProperty("FlatBackgroundCorrectionStop").value
        if (flat_background_correction_start == Property.EMPTY_DBL or
                        flat_background_correction_stop == Property.EMPTY_DBL):
            errors.update({"FlatBackgroundCorrectionStart": "The flat background correction start value needs "
                                                            "to be specified."})
            errors.update({"FlatBackgroundCorrectionStop": "The flat background correction stop value needs "
                                                           "to be specified."})

        if flat_background_correction_start > flat_background_correction_stop:
            errors.update({"FlatBackgroundCorrectionStart": "The flat background correction start value needs "
                                                            "to be smaller than the stop value."})
            errors.update({"FlatBackgroundCorrectionStop": "The flat background correction start value needs "
                                                           "to be smaller than the stop value."})

        # --------------------------
        # Rebin
        # --------------------------
        wavelength_low = self.getProperty("WavelengthLow").value
        wavelength_high = self.getProperty("WavelengthHigh").value
        wavelength_step = self.getProperty("WavelengthStep").value

        # The values need to be set
        if wavelength_low == Property.EMPTY_DBL:
            errors.update({"WavelengthLow": "The lower wavelength value needs to be specified."})
        if wavelength_high == Property.EMPTY_DBL:
            errors.update({"WavelengthHigh": "The upper wavelength value needs to be specified."})
        if wavelength_step == Property.EMPTY_DBL:
            errors.update({"WavelengthStep": "The wavelength step value needs to be specified."})

        if wavelength_low > wavelength_high:
            errors.update({"WavelengthLow": "The lower wavelength value needs to smaller than the"
                                            " upper wavelength value."})
            errors.update({"WavelengthHigh": "The lower wavelength value needs to smaller than the"
                                             " upper wavelength value."})
        return errors


# Register algorithm with Mantid
AlgorithmFactory.subscribe(SANSNormalizeToMonitor)
