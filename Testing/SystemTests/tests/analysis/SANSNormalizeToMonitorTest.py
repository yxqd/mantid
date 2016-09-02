import unittest
import stresstesting

import mantid
from SANS2.Common.SANSConstants import SANSConstants
from SANS2.Common.SANSFunctions import create_unmanaged_algorithm


def get_sample_monitor_workspace_for_normalize_to_monitor(use_peak=False):
    def set_up_data(data, use_prompt_peak):
        # From index 0 to 15  we set the data to 3. This is used as a flat background
        for index in range(16):
            data[index] = 3

        # From 16 to 100 set to 3. This simulates a measurement
        range_signal = range(16, 101)
        max_signal = len(range_signal) + 3
        for index in range(16, 101):
            data[index] = max_signal
            max_signal -= 1

        # Add a prompt peak
        if use_prompt_peak:
            for index in range(70, 80):
                data[index] = 100

    load_name = "LoadISISNexus"
    load_options = {"Filename": "LOQ74014",
                    SANSConstants.output_workspace: "normalize_to_monitor_output_workspace",
                    "LoadMonitors": "Separate"}
    load_alg = create_unmanaged_algorithm(load_name, **load_options)
    load_alg.execute()

    workspace = load_alg.getProperty(SANSConstants.output_monitor_workspace).value
    data_y0 = workspace.dataY(0)
    data_y1 = workspace.dataY(1)
    set_up_data(data_y0, use_peak)
    set_up_data(data_y1, use_peak)
    return workspace


class SANSNormalizeToMonitorTest(unittest.TestCase):
    monitor_workspace = None

    def _get_normalization_alg(self, workspace, incident_spectrum_number=None, scale_factor=None, prompt_start=None,
                               prompt_stop=None, flat_start=None, flat_stop=None, wavelength_low=None,
                               wavelength_high=None, wavelength_step=None, wavelength_step_type=None):
        normalize_name = "SANSNormalizeToMonitor"
        normalize_options = {SANSConstants.input_workspace: workspace,
                             SANSConstants.output_workspace: SANSConstants.dummy}
        if incident_spectrum_number:
            normalize_options.update({"IncidentMonitorSpectrumNumber": incident_spectrum_number})
        if scale_factor:
            normalize_options.update({"ScaleFactor": scale_factor})
        if prompt_start:
            normalize_options.update({"PromptPeakCorrectionStart": prompt_start})
        if prompt_stop:
            normalize_options.update({"PromptPeakCorrectionStop": prompt_stop})
        if flat_start:
            normalize_options.update({"FlatBackgroundCorrectionStart": flat_start})
        if flat_stop:
            normalize_options.update({"FlatBackgroundCorrectionStop": flat_stop})
        if wavelength_low:
            normalize_options.update({"WavelengthLow": wavelength_low})
        if wavelength_high:
            normalize_options.update({"WavelengthHigh": wavelength_high})
        if wavelength_step:
            normalize_options.update({"WavelengthStep": wavelength_step})
        if wavelength_step_type:
            normalize_options.update({"WavelengthStepType": wavelength_step_type})
        normalize_alg = create_unmanaged_algorithm(normalize_name, **normalize_options)
        return normalize_alg

    def test_that_invalid_prompt_peak_setting_raises(self):
        # Arrange
        if self.monitor_workspace is None:
            self.monitor_workspace = get_sample_monitor_workspace_for_normalize_to_monitor()

        prompt_peak_min = 7000
        prompt_peak_max = 1000
        normalize_alg = self._get_normalization_alg(workspace=self.monitor_workspace,
                                                    incident_spectrum_number=1, scale_factor=1.0,
                                                    prompt_start=prompt_peak_min, prompt_stop=prompt_peak_max,
                                                    flat_start=6000, flat_stop=16000,
                                                    wavelength_low=4., wavelength_high=16.,
                                                    wavelength_step=2., wavelength_step_type="LIN")

        # Act + Assert
        try:
            normalize_alg.execute()
            did_raise = False
        except RuntimeError:
            did_raise = True
        self.assertTrue(did_raise)

    def test_that_invalid_flat_background_time_setting_raises(self):
        # Arrange
        if self.monitor_workspace is None:
            self.monitor_workspace = get_sample_monitor_workspace_for_normalize_to_monitor()

        flat_start = 7000
        flat_stop = 1000
        normalize_alg = self._get_normalization_alg(workspace=self.monitor_workspace,
                                                    incident_spectrum_number=1, scale_factor=1.0,
                                                    prompt_start=1000, prompt_stop=2000,
                                                    flat_start=flat_start, flat_stop=flat_stop,
                                                    wavelength_low=4., wavelength_high=16.,
                                                    wavelength_step=2., wavelength_step_type="LIN")

        # Act + Assert
        try:
            normalize_alg.execute()
            did_raise = False
        except RuntimeError:
            did_raise = True
        self.assertTrue(did_raise)

    def test_that_invalid_wavelength_setting_raises(self):
        # Arrange
        if self.monitor_workspace is None:
            self.monitor_workspace = get_sample_monitor_workspace_for_normalize_to_monitor()

        wavelength_min = 8
        wavelength_max = 3
        normalize_alg = self._get_normalization_alg(workspace=self.monitor_workspace,
                                                    incident_spectrum_number=1, scale_factor=1.0,
                                                    flat_start=6000, flat_stop=12000,
                                                    wavelength_low=wavelength_min, wavelength_high=wavelength_max,
                                                    wavelength_step=2., wavelength_step_type="LIN")

        # Act + Assert
        try:
            normalize_alg.execute()
            did_raise = False
        except RuntimeError:
            did_raise = True
        self.assertTrue(did_raise)

    def test_that_workspace_with_prompt_peak_has_prompt_peak_removed(self):
        # Arrange
        if self.monitor_workspace is None:
            self.monitor_workspace = get_sample_monitor_workspace_for_normalize_to_monitor(True)

        normalize_alg = self._get_normalization_alg(workspace=self.monitor_workspace,
                                                    incident_spectrum_number=1, scale_factor=1.0,
                                                    prompt_start=19900, prompt_stop=25000,
                                                    flat_start=3500, flat_stop=5100,
                                                    wavelength_low=4., wavelength_high=14.5,
                                                    wavelength_step=.3, wavelength_step_type="LIN")
        # Act
        normalize_alg.execute()

        # Assert
        output_workspace = normalize_alg.getProperty(SANSConstants.output_workspace).value
        data_y = output_workspace.dataY(0)
        # Removal of our artificial prompt peak will lead to a monotonously decreasing data.
        for index in range(0, len(data_y) - 1):
            index_advanced = index + 1
            self.assertTrue(data_y[index] >= data_y[index_advanced])


class SANSNormalizeToMonitorRunnerTest(stresstesting.MantidStressTest):
    def __init__(self):
        stresstesting.MantidStressTest.__init__(self)
        self._success = False

    def runTest(self):
        suite = unittest.TestSuite()
        suite.addTest(unittest.makeSuite(SANSNormalizeToMonitorTest, 'test'))
        runner = unittest.TextTestRunner()
        res = runner.run(suite)
        if res.wasSuccessful():
            self._success = True

    def requiredMemoryMB(self):
        return 2000

    def validate(self):
        return self._success


if __name__ == '__main__':
    unittest.main()
