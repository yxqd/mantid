# pylint: disable=invalid-name

""" SANSCreateWavelengthAndPixelAdjustment algorithm creates workspaces for pixel adjustment
    and wavelength adjustment.
"""

from mantid.kernel import (Direction, IntBoundedValidator, FloatBoundedValidator, StringListValidator, Property,
                           PropertyManagerProperty)
from mantid.api import (DataProcessorAlgorithm, MatrixWorkspaceProperty, AlgorithmFactory, PropertyMode,
                        FileProperty, FileAction, Progress)
from SANS2.State.SANSStateBase import create_deserialized_sans_state_from_property_manager
from SANS2.Common.SANSEnumerations import (RangeStepType)
from SANS2.Common.SANSConstants import SANSConstants
from SANS2.Common.SANSFunctions import create_unmanaged_algorithm


class SANSCreateWavelengthAndPixelAdjustment(DataProcessorAlgorithm):
    def category(self):
        return 'SANS\\Adjust'

    def summary(self):
        return 'Calculates wavelength adjustment and pixel adjustment workspaces.'

    def PyInit(self):
        # State
        self.declareProperty(PropertyManagerProperty('SANSState'),
                             doc='A property manager which fulfills the SANSState contract.')
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

        # The component
        self.declareProperty('Component', '', direction=Direction.Input, doc='Component which is being investigated.')

    def PyExec(self):
        # Read the state
        state_property_manager = self.getProperty("SANSState").value
        state = create_deserialized_sans_state_from_property_manager(state_property_manager)
        wavelength_and_pixel_adjustment_state = state.adjustment.wavelength_and_pixel_adjustment

        # Get the wavelength adjustment workspace
        transmission_workspace = self.getProperty("Transmission").value
        monitor_normalization_workspace = self.getProperty("MonitorNormalization").value

        wavelength_adjustment_file = wavelength_and_pixel_adjustment_state.wavelength_adjustment_file
        rebin_string = self._get_rebin_string(wavelength_and_pixel_adjustment_state)
        wavelength_adjustment_workspace = self._get_wavelength_adjustment_workspace(wavelength_adjustment_file,
                                                                                    transmission_workspace,
                                                                                    monitor_normalization_workspace,
                                                                                    rebin_string)

        # Get the pixel adjustment workspace
        component = self.getProperty("Component").value
        pixel_adjustment_file = wavelength_and_pixel_adjustment_state.pixel_adjustment_file
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

    def _get_rebin_string(self, wavelength_and_pixel_adjustment_state):
        wavelength_low = wavelength_and_pixel_adjustment_state.wavelength_low
        wavelength_high = wavelength_and_pixel_adjustment_state.wavelength_high
        wavelength_step = wavelength_and_pixel_adjustment_state.wavelength_step
        wavelength_step_type = -1.0 if wavelength_and_pixel_adjustment_state.wavelength_step_type \
                                       is RangeStepType.Log else 1.0

        # Create a rebin string from the wavelength information
        wavelength_step *= wavelength_step_type
        return str(wavelength_low) + "," + str(wavelength_step) + "," + str(wavelength_high)

    def validateInputs(self):
        errors = dict()
        # Check that the input can be converted into the right state object
        state_property_manager = self.getProperty("SANSState").value
        try:
            state = create_deserialized_sans_state_from_property_manager(state_property_manager)
            state.property_manager = state_property_manager
            state.validate()
        except ValueError as err:
            errors.update({"SANSSMove": str(err)})
        return errors

# Register algorithm with Mantid
AlgorithmFactory.subscribe(SANSCreateWavelengthAndPixelAdjustment)
