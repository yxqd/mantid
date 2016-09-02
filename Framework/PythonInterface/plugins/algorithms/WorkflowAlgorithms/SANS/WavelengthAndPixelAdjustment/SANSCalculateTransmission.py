# pylint: disable=invalid-name

""" SANSCalculateTransmission algorithm calculates the transmission correction of a SANS workspace."""

from mantid.kernel import (Direction, StringArrayProperty, StringListValidator, Property,
                           PropertyManagerProperty, FloatBoundedValidator)
from mantid.api import (DataProcessorAlgorithm, MatrixWorkspaceProperty, AlgorithmFactory, PropertyMode,
                        FileProperty, FileAction, Progress)
from SANS2.Common.SANSConstants import SANSConstants
from SANS2.Common.SANSFunctions import create_unmanaged_algorithm
from SANS2.WavelengthAndPixelAdjustment.CalculateTransmissionHelper import (get_detector_id_for_spectrum_number,
                                                                            get_workspace_indices_for_monitors,
                                                                            apply_flat_background_correction_to_monitors,
                                                                            apply_flat_background_correction_to_detectors,
                                                                            get_region_of_interest)


class SANSCalculateTransmission(DataProcessorAlgorithm):
    def category(self):
        return 'SANS\\Transmission'

    def summary(self):
        return 'Calculates the transmission for a SANS reduction.'

    def PyInit(self):
        # ---------------------------
        # Workspaces
        # ---------------------------
        # Input workspace in TOF
        self.declareProperty(MatrixWorkspaceProperty("TransmissionWorkspace", '',
                                                     optional=PropertyMode.Mandatory, direction=Direction.Input),
                             doc='The transmission workspace in time-of-flight units.')
        self.declareProperty(MatrixWorkspaceProperty("DirectWorkspace", '',
                                                     optional=PropertyMode.Mandatory, direction=Direction.Input),
                             doc='The direct workspace in time-of-flight units.')

        # Output workspace
        self.declareProperty(MatrixWorkspaceProperty(SANSConstants.output_workspace, '', direction=Direction.Output),
                             doc='A workspace in units of wavelength.')

        # ---------------------------
        # The incident monitor
        # ---------------------------
        self.declareProperty('IncidentMonitorSpectrumNumber', defaultValue=Property.EMPTY_INT,
                             direction=Direction.Input,
                             doc='The spectrum number of the incident monitor')

        # ---------------------------
        # Transmission selection
        # ---------------------------
        self.declareProperty('TransmissionMonitor', defaultValue=Property.EMPTY_INT,
                             direction=Direction.Input,
                             doc='The spectrum number of the transmission monitor')
        self.declareProperty('TransmissionRadius', defaultValue=Property.EMPTY_DBL,
                             direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='The workspace index of the transmission monitor')
        self.declareProperty(StringArrayProperty("TransmissionROIFiles",
                                                 direction=Direction.Input),
                             doc="Comma separated list transmission region of interest files.")
        self.declareProperty(StringArrayProperty("TransmissionMaskFiles",
                                                 direction=Direction.Input),
                             doc="Comma separated list transmission mask files.")
        self.setPropertyGroup("TransmissionMonitor", 'Transmission Settings')
        self.setPropertyGroup("TransmissionRadius", 'Transmission Settings')
        self.setPropertyGroup("TransmissionROIFiles", 'Transmission Settings')
        self.setPropertyGroup("TransmissionMaskFiles", 'Transmission Settings')

        # ---------------------------
        # Prompt Peak correction.
        # ---------------------------
        self.declareProperty('PromptPeakCorrectionStart', defaultValue=Property.EMPTY_DBL,
                             direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='The start time of the prompt peak.')
        self.declareProperty('PromptPeakCorrectionStop', defaultValue=Property.EMPTY_DBL,
                             direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='The stop time of the prompt peak.')
        self.setPropertyGroup("PromptPeakCorrectionStart", 'Prompt Peak Correction')
        self.setPropertyGroup("PromptPeakCorrectionStop", 'Prompt Peak Correction')

        # ---------------------------
        # Flat background settings
        # ----------------------------
        self.declareProperty('FlatBackgroundCorrectionROIStart', defaultValue=Property.EMPTY_DBL,
                             validator=FloatBoundedValidator(0.0),
                             direction=Direction.Input,
                             doc='The start time of the flat background correction if a region of interest is used.')
        self.declareProperty('FlatBackgroundCorrectionROIStop', defaultValue=Property.EMPTY_DBL,
                             direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='The stop time of the flat background correction if a region of interest is used.')
        self.declareProperty(PropertyManagerProperty('FlatBackgroundCorrectionMonitors'),
                             doc='A property manager which contains the spectrum number (as a string) of monitors '
                                 'as a key and an array with exactly two entries as the value. The two entries are the '
                                 'start and stop time of the flat background correction of that particular monitor')
        self.setPropertyGroup("FlatBackgroundCorrectionROIStart", 'Flat Background Correction')
        self.setPropertyGroup("FlatBackgroundCorrectionROIStop", 'Flat Background Correction')
        self.setPropertyGroup("FlatBackgroundCorrectionMonitors", 'Flat Background Correction')

        # ---------------------------
        # Wavelength conversion
        # ---------------------------
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
        self.setPropertyGroup("WavelengthLow", 'Convert to Wavelength Settings')
        self.setPropertyGroup("WavelengthHigh", 'Convert to Wavelength Settings')
        self.setPropertyGroup("WavelengthStep", 'Convert to Wavelength Settings')
        self.setPropertyGroup("WavelengthStepType", 'Convert to Wavelength Settings')
        self.setPropertyGroup("RebinMode", 'Convert to Wavelength Settings')

        # ---------------------------
        # Fitting
        # ---------------------------
        allowed_fit_method = StringListValidator(["Log", "Linear", "Polynomial"])
        self.declareProperty('FitMethod', "LIN", validator=allowed_fit_method, direction=Direction.Input,
                             doc='The fit method for the transmission calculation.')
        self.declareProperty('PolynomialOrder', defaultValue=Property.EMPTY_INT, direction=Direction.Input,
                             doc='The polynomial order of a polynomial fit.  It is considered only '
                                 'for FitMethod=Polynomial')
        # The output will be returned in wavelength units. These properties are for the wavelength rebinning.
        self.declareProperty('FitWavelengthLow', defaultValue=Property.EMPTY_DBL, direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='The low value of the wavelength range for the fit.')
        self.declareProperty('FitWavelengthHigh', defaultValue=Property.EMPTY_DBL, direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='The high value of the wavelength binning for the fit.')
        self.declareProperty('FitWavelengthStep', defaultValue=Property.EMPTY_DBL, direction=Direction.Input,
                             validator=FloatBoundedValidator(0.0),
                             doc='The step size of the wavelength binning for the fit.')
        allowed_step_types = StringListValidator(["LOG", "LIN"])
        self.declareProperty('FitWavelengthStepType', "LIN", validator=allowed_step_types, direction=Direction.Input,
                             doc='The step type for rebinning for the fit.')
        self.setPropertyGroup("FitMethod", 'Fit')
        self.setPropertyGroup("PolynomialOrder", 'Fit')
        self.setPropertyGroup("FitWavelengthLow", 'Fit')
        self.setPropertyGroup("FitWavelengthHigh", 'Fit')
        self.setPropertyGroup("FitWavelengthStep", 'Fit')
        self.setPropertyGroup("FitWavelengthStepType", 'Fit')

    def PyExec(self):
        # The calculation of the transmission has the following steps:
        # 1. Get all spectrum numbers which take part in the transmission calculation
        # 2. Clean up the transmission and direct workspaces, ie peak prompt correction, flat background calculation,
        #    wavelength conversion and rebinning of the data.
        # 3. Run the CalculateTransmission algorithm
        transmission_workspace = self.getProperty("TransmissionWorkspace").value
        direct_workspace = self.getProperty("DirectWorkspace").value
        incident_monitor_spectrum_number = self.getProperty("IncidentMonitorSpectrumNumber").value

        # 1. Get relevant spectra
        detector_id_incident_monitor = get_detector_id_for_spectrum_number(transmission_workspace,
                                                                           incident_monitor_spectrum_number)
        detector_ids_roi, detector_id_transmission_monitor = self._get_detector_ids_for_transmission_calculation(
                                                                  transmission_workspace)
        all_detector_ids = [detector_id_incident_monitor]
        if detector_ids_roi:
            all_detector_ids.extend(detector_ids_roi)
        elif detector_id_transmission_monitor:
            all_detector_ids.append(detector_id_transmission_monitor)
        else:
            raise RuntimeError("SANSCalculateTransmission: No region of interest or transmission monitor selected.")

        # 2. Clean transmission data
        transmission_workspace = self._get_corrected_wavelength_workspace(transmission_workspace, all_detector_ids)
        direct_workspace = self._get_corrected_wavelength_workspace(direct_workspace, all_detector_ids)

        # 3. Fit
        fitted_transmission_workspace, unfitted_transmission_workspace = \
            self._perform_fit(transmission_workspace, direct_workspace, detector_ids_roi,
                              detector_id_transmission_monitor, detector_id_incident_monitor)

        # Set the output
        self.setProperty(SANSConstants.output_workspace, fitted_transmission_workspace)
        self.setProperty("UnfittedData", unfitted_transmission_workspace)

    def _perform_fit(self, transmission_workspace, direct_workspace,
                     transmission_roi_detector_ids, transmission_monitor_detector_id, incident_monitor_detector_id):
        wavelength_low = self.getProperty("FitWavelengthLow").value
        wavelength_high = self.getProperty("FitWavelengthHigh").value
        wavelength_step = self.getProperty("FitWavelengthStep").value
        wavelength_step_type = self.getProperty("FitWavelengthStepType").value
        prefix = 1.0 if wavelength_step_type == "LIN" else -1.0
        wavelength_step *= prefix
        rebin_params = str(wavelength_low) + "," + str(wavelength_step) + "," + str(wavelength_high)

        trans_name = "CalculateTransmission"
        trans_options = {"SampleRunWorkspace": transmission_workspace,
                         "DirectRunWorkspace": direct_workspace,
                         SANSConstants.output_workspace: SANSConstants.dummy,
                         "IncidentBeamMonitor": incident_monitor_detector_id,
                         "RebinParams": rebin_params,
                         "OutputUnfittedData": True}

        # If we have a region of interest we use it else we use the transmission monitor
        if transmission_roi_detector_ids:
            trans_options.update({"TransmissionROI": transmission_roi_detector_ids})
        else:
            trans_options.update({"TransmissionMonitor", transmission_monitor_detector_id})

        fit_method = self.getProperty("FitMethod").value
        trans_options.update({"FitMethod": fit_method})
        if fit_method == "Polynomial":
            polynomial_order = self.getProperty("PolynomialOrder").value
            trans_options.update({"PolynomialOrder": polynomial_order})

        trans_alg = create_unmanaged_algorithm(trans_name, **trans_options)
        trans_alg.execute()

        fitted_transmission_workspace = self.getProperty(SANSConstants.output_workspace).value
        unfitted_transmission_workspace = self.getProperty("UnfittedData").value

        # Set the y label correctly for the fitted and unfitted transmission workspaces
        y_unit_label_transmission_ratio = "Transmission"
        if fitted_transmission_workspace:
            fitted_transmission_workspace.setYUnitLabel(y_unit_label_transmission_ratio)
        if unfitted_transmission_workspace:
            unfitted_transmission_workspace.setYUnitLabel(y_unit_label_transmission_ratio)
        return fitted_transmission_workspace, unfitted_transmission_workspace

    def _get_detector_ids_for_transmission_calculation(self, transmission_workspace):
        """
        Get the detector ids which participate in the transmission calculation.

        This can come either from a ROI/MASK/RADIUS selection or from a transmission monitor, not both.
        :param transmission_workspace: the transmission workspace.
        :return: a list of detector ids for ROI and a detector id for the transmission monitor, either can be None
        """
        transmission_radius = self.geProperty("TransmissionRadius").value
        if transmission_radius == Property.EMPTY_DBL:
            transmission_radius = None
        transmission_roi = self.getProperty("TransmissionROIFiles").value
        transmission_mask = self.getProperty("TransmissionMaskFiles").value
        detector_ids_roi = get_region_of_interest(transmission_workspace, transmission_radius, transmission_roi,
                                                  transmission_mask)

        transmission_monitor_spectrum_number = self.getProperty("TransmissionMonitor").value
        detector_id_transmission_monitor = get_detector_id_for_spectrum_number(transmission_monitor_spectrum_number)

        return detector_ids_roi, detector_id_transmission_monitor

    def _get_corrected_wavelength_workspace(self, workspace, detector_ids):
        """
        Performs a prompt peak correction, a background correction, converts to wavelength and rebins.

        :param workspace: the workspace which is being corrected.
        :param detector_ids: a list of relevant detector ids
        :return:  a corrected workspace.
        """
        # Extract the relevant spectra. These include
        # 1. The incident monitor spectrum
        # 2. The transmission spectra, be it monitor or ROI based.
        # A previous implementation of this code had a comment which suggested
        # that we have to exclude unused spectra as the interpolation runs into
        # problems if we don't.
        extract_name = "ExtractSpectra"
        extract_options = {SANSConstants.input_workspace: workspace,
                           SANSConstants.output_workspace: SANSConstants.dummy,
                           "DetectorList": detector_ids}
        extract_alg = create_unmanaged_algorithm(extract_name, **extract_options)
        extract_alg.execute()
        workspace = extract_alg.getProperty(SANSConstants.output_workspace).value

        # Make sure that we still have spectra in the workspace
        if workspace.getNumberHistograms() == 0:
            raise RuntimeError("SANSCalculateTransmissionCorrection: The transmission workspace does "
                               "not seem to have any spectra.")

        # ----------------------------------
        # Perform the prompt peak correction
        # ----------------------------------
        prompt_peak_correction_start = self.getProperty("PromptPeakCorrectionStart").value
        prompt_peak_correction_stop = self.getProperty("PromptPeakCorrectionStop").value
        workspace = self._perform_prompt_peak_correction(workspace, prompt_peak_correction_start,
                                                         prompt_peak_correction_stop)

        # ---------------------------------------
        # Perform the flat background correction
        # ---------------------------------------
        # The flat background correction has two parts:
        # 1. Corrections on monitors
        # 2. Corrections on regular detectors

        # Monitor flat background correction
        workspace_indices_of_monitors = list(get_workspace_indices_for_monitors(workspace))
        flat_background_monitors = self.getProperty("FlatBackgroundCorrectionMonitors").value
        workspace = apply_flat_background_correction_to_monitors(workspace,  workspace_indices_of_monitors,
                                                                 flat_background_monitors)

        # Detector flat background correction
        flat_background_correction_start = self.getProperty("FlatBackgroundCorrectionROIStart").value
        flat_background_correction_stop = self.getProperty("FlatBackgroundCorrectionROIStop").value
        workspace = apply_flat_background_correction_to_detectors(workspace, flat_background_correction_start,
                                                                  flat_background_correction_stop)

        # ---------------------------------------
        # Convert to wavelength and rebin
        # ---------------------------------------
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

    def _perform_prompt_peak_correction(self, workspace, prompt_peak_correction_start, prompt_peak_correction_stop):
        """
        Prompt peak correction is performed if it is explicitly set by the user.

        :param workspace: the workspace to correct.
        :param prompt_peak_correction_start: the start time for the prompt peak correction.
        :param prompt_peak_correction_stop: the stop time for the prompt peak correction.
        :return: a corrected workspace.
        """
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
            workspace = remove_alg.getProperty(SANSConstants.output_workspace)
        return workspace

    def validateInputs(self):
        errors = dict()

        transmission_workspace = self.getProperty("TransmissionWorkspace").value

        # ------------------------------
        # IncidentMonitorSpectrumNumber
        # ------------------------------
        incident_transmission_spectrum = self.getProperty("IncidentMonitorSpectrumNumber").value
        # Make sure that the spectrum number exists on the transmission workspace
        try:
            transmission_workspace.getIndexFromSpectrumNumber(incident_transmission_spectrum)
        except RuntimeError:
            errors.update({"IncidentMonitorSpectrumNumber": "The spectrum number for the incident monitor spectrum "
                                                            "does not seem to exist for the transmission workspace."})

        # ------------------------------
        # Transmission Selection
        # ------------------------------
        trans_monitor_spectrum_number = self.getProperty("TransmissionMonitor").value
        trans_radius = self.getProperty("TransmissionRadius").value
        trans_roi = self.getProperty("TransmissionROIFiles").value
        trans_mask = self.getProperty("TransmissionMaskFiles").value

        # Either the transmission monitor or a ROI selection must be made
        roi_selection = trans_radius is None and trans_roi is None and trans_mask is None
        if trans_monitor_spectrum_number is None and roi_selection is None:
            errors.update({"TransmissionMonitor": "Either a transmission monitor or a ROI "
                                                  "selection must be performed."})
            errors.update({"TransmissionRadius": "Either a transmission monitor or a ROI "
                                                 "selection must be performed."})
            errors.update({"TransmissionROIFiles": "Either a transmission monitor or a ROI "
                                                   "selection must be performed."})
            errors.update({"TransmissionMaskFiles": "Either a transmission monitor or a ROI "
                                                    "selection must be performed."})

        # ---------------------------
        # Prompt Peak correction.
        # ---------------------------
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

        # ---------------------------
        # Flat background settings
        # ----------------------------
        flat_background_correction_start = self.getProperty("FlatBackgroundCorrectionROIStart").value
        flat_background_correction_stop = self.getProperty("FlatBackgroundCorrectionROIStop").value
        if (flat_background_correction_start == Property.EMPTY_DBL or
                        flat_background_correction_stop == Property.EMPTY_DBL):
            errors.update({"FlatBackgroundCorrectionROIStart": "The flat background correction start value needs "
                                                               "to be specified."})
            errors.update({"FlatBackgroundCorrectionROIStop": "The flat background correction stop value needs "
                                                              "to be specified."})

        if flat_background_correction_start > flat_background_correction_stop:
            errors.update({"FlatBackgroundCorrectionROIStart": "The flat background correction start value needs "
                                                               "to be smaller than the stop value."})
            errors.update({"FlatBackgroundCorrectionROIStop": "The flat background correction start value needs "
                                                              "to be smaller than the stop value."})

        # Ensure that the monitor flat_background has only key values which are convertible to integers and exactly two
        # entries where the first entry is smaller than the second
        flat_background_correction = self.getProperty("FlatBackgroundCorrectionMonitors").value
        keys = flat_background_correction.keys()
        for key in keys:
            try:
                int(key)
            except ValueError:
                errors.update({"FlatBackgroundCorrectionMonitors": "The key {0} cannot be converted to an "
                                                                   "integer".format(key)})
            value = flat_background_correction[key]

            try:
                # Check that it has two entries
                if len(value) != 2:
                    errors.update({"FlatBackgroundCorrectionMonitors": "The list {0} needs to have exactly two "
                                                                       "entries.".format(key)})
                # Check that the first entry is smaller than the second entry
                if value[0] > value[1]:
                    errors.update({"FlatBackgroundCorrectionMonitors": "The first entry {0} is larger than the second "
                                                                       "entry {1}, but it needs to be "
                                                                       "smaller.".format(value[0], value[1])})
            except TypeError:
                flat_background_correction[key] = errors.update({"FlatBackgroundCorrectionMonitors": "The value must "
                                                            "be a collection, but it is a {0}.".format(type(value))})

        # ---------------------------
        # Wavelength conversion
        # ---------------------------
        wavelength_low = self.getProperty("WavelengthLow").value
        wavelength_high = self.getProperty("WavelengthHigh").value
        # Make sure that the lower wavelength is larger than the upper wavelength
        if wavelength_low > wavelength_high:
            errors.update({"WavelengthLow": "The lower wavelength setting is larger than the higher "
                                            "wavelength setting."})
            errors.update({"WavelengthHigh": "The lower wavelength setting is larger than the higher "
                                             "wavelength setting."})

        # ---------------------------
        # Fit
        # ---------------------------
        fit_wavelength_low = self.getProperty("FitWavelengthLow").value
        fit_wavelength_high = self.getProperty("FitWavelengthHigh").value
        # Make sure that the lower wavelength is larger than the upper wavelength
        if fit_wavelength_low > fit_wavelength_high:
            errors.update({"FitWavelengthLow": "The lower wavelength setting is larger than the higher "
                                               "wavelength setting."})
            errors.update({"FitWavelengthHigh": "The lower wavelength setting is larger than the higher "
                                                "wavelength setting."})
        return errors


# Register algorithm with Mantid
AlgorithmFactory.subscribe(SANSCalculateTransmission)
