# pylint: disable=invalid-name

""" SANSCreateWavelengthAndPixelAdjustment algorithm creates workspaces for pixel adjustment
    and wavelength adjustment.
"""

from mantid.kernel import (Direction, IntBoundedValidator, FloatBoundedValidator, StringListValidator, Property)
from mantid.api import (DataProcessorAlgorithm, MatrixWorkspaceProperty, AlgorithmFactory, PropertyMode,
                        FileProperty, FileAction, Progress)
from mantid.dataobjects import EventWorkspace

from SANS2.Common.SANSConstants import SANSConstants
from SANS2.Common.SANSFunctions import create_unmanaged_algorithm


class SANSCreateWavelengthAndPixelAdjustment(DataProcessorAlgorithm):
    def category(self):
        return 'SANS\\Normalize'

    def summary(self):
        return 'Calculates wavelength adjustment and pixel adjustment workspaces.'

    def PyInit(self):
        # Input workspaces
        self.declareProperty(MatrixWorkspaceProperty("Transmission", '',
                                                     optional=PropertyMode.Mandatory, direction=Direction.Input),
                             doc='The calculated transmission workspace in wavelength units.')
        self.declareProperty(MatrixWorkspaceProperty("MonitorNormalization", '',
                                                     optional=PropertyMode.Mandatory, direction=Direction.Input),
                             doc='The monitor normalization workspace in wavelength units.')

        # Output workspace
        self.declareProperty(MatrixWorkspaceProperty("OutputWorkspaceWavelengthAdjustment", '',
                                                     direction=Direction.Output),
                             doc='A wavelength adjustment output workspace.')
        self.declareProperty(MatrixWorkspaceProperty("OutputWorkspacePixelAdjustment", '',
                                                     direction=Direction.Output),
                             doc='A pixel adjustment output workspace.')

        # The adjustment files for wavelength adjustments and pixel adjustments
        self.declareProperty(FileProperty("WavelengthAdjustmentFile", "", action=FileAction.Load, extensions=[".txt"]))
        self.declareProperty(FileProperty("PixelAdjustmentFile", "", action=FileAction.Load, extensions=[".txt"]))

        # The component
        self.declareProperty('Component', '', direction=Direction.Input, doc='Component which is being investigated.')

        # Wavelength rebin values
        self.declareProperty('WavelengthLow', defaultValue=Property.EMPTY_DBL, direction=Direction.Input,
                             doc='The low value of the wavelength binning.')
        self.declareProperty('WavelengthHigh', defaultValue=Property.EMPTY_DBL, direction=Direction.Input,
                             doc='The high value of the wavelength binning.')
        self.declareProperty('WavelengthStep', defaultValue=Property.EMPTY_DBL, direction=Direction.Input,
                             doc='The step size of the wavelength binning.')
        allowed_step_types = StringListValidator(["LOG", "LIN"])
        self.declareProperty('WavelengthStepType', "LIN", validator=allowed_step_types, direction=Direction.Input,
                             doc='The step type for rebinning.')

    def PyExec(self):
        # Get the wavelength adjustment workspace
        wavelength_adjustment_file = self.getProperty("WavelengthAdjustmentFile").value
        transmission_workspace = self.getProperty("Transmission").value
        monitor_normalization_workspace = self.getProperty("MonitorNormalization").value

        rebin_string = self._get_rebin_string()
        wavelength_adjustment_workspace = self._get_wavelength_adjustment_workspace(wavelength_adjustment_file,
                                                                                    transmission_workspace,
                                                                                    monitor_normalization_workspace,
                                                                                    rebin_string)

        # Get the pixel adjustment workspace
        component = self.getProperty("Component").value
        pixel_adjustment_file = self.getProperty("PixelAdjustmentFile").value
        pixel_adjustment_workspace = self._get_pixel_adjustment_workspace(pixel_adjustment_file, component)

        # Set the output
        self.setProperty("OutputWorkspaceWavelengthAdjustment", wavelength_adjustment_workspace)
        self.setProperty("OutputWorkspacePixelAdjustment", pixel_adjustment_workspace)

    def _get_wavelength_adjustment_workspace(self, wavelength_adjustment_file, transmission_workspace,
                                             monitor_normalization_workspace, rebin_string):
        """
        This creates a workspace with wavelength adjustments, ie this will be a correction for the bins, but it will
        be the same for all pixels. This is essentially the product of several workspaces.
        The participating workspaces are:
        1. A workspace loaded from a calibration file
        2. The workspace resulting from the transmission calcuation (using SANSCalculateTransmission)
        3. The workspace resulting from the monitor normalization

        :param wavelength_adjustment_file: the file path to the wavelength adjustment file
        :param transmission_workspace: the calculated transmission workspace
        :param monitor_normalization_workspace: the monitor normalization workspace
        :param rebin_string: the parameters for rebinning
        :return: a general wavelength adjustment workspace
        """
        # Get the wavelength correction workspace from the file
        wavelength_correction_workspace_from_file = self._load_wavelength_correction_file(wavelength_adjustment_file)

        wavelength_adjustment_workspaces = [wavelength_correction_workspace_from_file, monitor_normalization_workspace,
                                            transmission_workspace]

        # Multiply all workspaces
        wavelength_adjustment_workspace = None
        for workspace in wavelength_adjustment_workspaces:
            if workspace is not None:
                # First we need to change the binning such that is matches the binning of the main data workspace
                rebin_name = "Rebin"
                rebin_options = {SANSConstants.input_workspace: workspace,
                                 "Params": rebin_string,
                                 SANSConstants.output_workspace: SANSConstants.dummy}
                rebin_alg = create_unmanaged_algorithm(rebin_name, **rebin_options)
                rebin_alg.execute()
                rebinned_workspace = rebin_alg.getProperty(SANSConstants.output_workspace).value

                if wavelength_adjustment_workspace is None:
                    wavelength_adjustment_workspace = rebinned_workspace
                else:
                    multiply_name = "Multiply"
                    multiply_options = {"LHSWorkspace": rebinned_workspace,
                                        "RHSWorkspace": wavelength_adjustment_workspace,
                                        SANSConstants.output_workspace: SANSConstants.dummy}
                    multiply_alg = create_unmanaged_algorithm(multiply_name, **multiply_options)
                    multiply_alg.execute()
                    wavelength_adjustment_workspace = multiply_alg.getProperty(SANSConstants.output_workspace).value
            return wavelength_adjustment_workspace

    def _load_wavelength_correction_file(self, file_name):
        correction_workspace = None
        if file_name:
            load_name = "LoadRKH"
            load_option = {"Filename": file_name,
                           SANSConstants.output_workspace: SANSConstants.dummy,
                            "FirstColumnValue": "Wavelength"}
            load_alg = create_unmanaged_algorithm(load_name, **load_option)
            load_alg.execute()
            output_workspace = load_alg.getProperty(SANSConstants.output_workspace).value
            # We require HistogramData and not PointData
            if not output_workspace.isHistogramData():
                convert_name = "ConvertToHistogram"
                convert_options = {SANSConstants.input_workspace: output_workspace,
                                   SANSConstants.output_workspace: SANSConstants.dummy}
                convert_alg = create_unmanaged_algorithm(convert_name, **convert_options)
                convert_alg.execute()
                output_workspace = convert_alg.getProperty(SANSConstants.output_workspace).value
            correction_workspace = output_workspace
        return correction_workspace

    def _get_pixel_adjustment_workspace(self, pixel_adjustment_file, component):
        """
        This get the pixel-by-pixel adjustment of the workspace

        :param pixel_adjustment_file: full file path to the pixel adjustment file
        :param component: the component which is currently being investigated
        :return: the pixel adjustment workspace
        """

        if pixel_adjustment_file:
            load_name = "LoadRKH"
            load_options = {SANSConstants.file_name: pixel_adjustment_file,
                            SANSConstants.output_workspace: SANSConstants.dummy,
                            "FirstColumnValue": "SpectrumNumber"}
            load_alg = create_unmanaged_algorithm(load_name, **load_options)
            load_alg.execute()
            output_workspace = load_alg.getProperty(load_name, **load_options).value

            # TODO: Crop to the detector under investigation
            crop_name = "SANSCrop"
            crop_options = {SANSConstants.input_workspace: output_workspace,
                            SANSConstants.output_workspace: SANSConstants.dummy,
                            "Component": component}
            crop_alg = create_unmanaged_algorithm(crop_name, **crop_options)
            crop_alg.execute()
            pixel_adjustment_workspace = crop_alg.getProperty(SANSConstants.output_workspace).value
        else:
            pixel_adjustment_workspace = None
        return pixel_adjustment_workspace

    def _get_rebin_string(self):
        wavelength_low = self.getProperty("WavelengthLow").value
        wavelength_high = self.getProperty("WavelengthHigh").value
        wavelength_step = self.getProperty("WavelengthStep").value
        wavelength_step_type = self.getProperty("WavelengthStepType").value

        # Create a rebin string from the wavelength information
        wavelength_step = -1*wavelength_step if wavelength_step_type == "LOG" else wavelength_step
        return str(wavelength_low) + "," + str(wavelength_step) + "," + str(wavelength_high)

    def validateInputs(self):
        errors = dict()
        # Check the wavelength
        wavelength_low = self.getProperty("WavelengthLow").value
        wavelength_high = self.getProperty("WavelengthHigh").value
        if wavelength_low is not None and wavelength_high is not None and wavelength_low > wavelength_high:
            errors.update({"WavelengthLow": "The lower wavelength setting needs to be smaller "
                                            "than the higher wavelength setting."})

        if wavelength_low is not None and wavelength_low < 0:
            errors.update({"WavelengthLow": "The wavelength cannot be smaller than 0."})

        if wavelength_high is not None and wavelength_high < 0:
            errors.update({"WavelengthHigh": "The wavelength cannot be smaller than 0."})

        wavelength_step = self.getProperty("WavelengthStep").value
        if wavelength_step is not None and wavelength_step < 0:
            errors.update({"WavelengthStep": "The wavelength step cannot be smaller than 0."})

        # Check the workspace
        workspace = self.getProperty(SANSConstants.input_workspace).value
        rebin_mode = self.getProperty("RebinMode").value
        if rebin_mode == "InterpolatingRebin" and isinstance(workspace, EventWorkspace):
            errors.update({"RebinMode": "An interpolating rebin cannot be applied to an EventWorkspace."})
        return errors

# Register algorithm with Mantid
AlgorithmFactory.subscribe(SANSCreateWavelengthAndPixelAdjustment)
