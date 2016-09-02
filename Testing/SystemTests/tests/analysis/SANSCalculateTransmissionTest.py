import unittest
import stresstesting

import mantid
from SANS2.Common.SANSConstants import SANSConstants
from SANS2.Common.SANSFunctions import create_unmanaged_algorithm

class SANSCalculateTransmissionTest(unittest.TestCase):
    def _perform_calculation(self, trans_workspace, direct_workspace, incid_mon=None, trans_mon=None, trans_radius=None,
                             trans_roi=None, trans_mask=None, prompt_start=None, prompt_stop=None, flat_start=None,
                             flat_stop=None, flat_mon=None, wav_low=None, wav_high=None, wav_step=None,
                             wav_step_type=None, rebin_mode=None, fit_method=None, poly_order=None, fit_wav_low=None,
                             fit_wav_high=None, fit_wav_step=None, fit_wav_step_type=None):
        transmission_name = "SANSCalculateTransmission"
        transmission_options = {"TransmissionWorkspace": trans_workspace,
                                "DirectWorkspace": direct_workspace,
                                SANSConstants.output_workspace: SANSConstants.dummy}
        if incid_mon:
            transmission_options.update({"IncidentMonitorSpectrumNumber": incid_mon})
        # --------------
        if trans_mon:
            transmission_options.update({"TransmissionMonitor": trans_mon})
        if trans_radius:
            transmission_options.update({"TransmissionRadius": trans_radius})
        if trans_roi:
            transmission_options.update({"TransmissionROIFiles": trans_roi})
        if trans_mask:
            transmission_options.update({"TransmissionMaskFiles": trans_mask})
        # --------------
        if prompt_start:
            transmission_options.update({"PromptPeakCorrectionStart": prompt_start})
        if prompt_stop:
            transmission_options.update({"PromptPeakCorrectionStop": prompt_stop})
        # --------------
        if flat_start:
            transmission_options.update({"FlatBackgroundCorrectionROIStart": flat_start})
        if flat_stop:
            transmission_options.update({"PromptPeakCorrectionStop": flat_stop})
        if flat_mon:
            transmission_options.update({"FlatBackgroundCorrectionMonitors": flat_mon})
        # --------------
        if wav_low:
            transmission_options.update({"WavelengthLow": wav_low})
        if wav_high:
            transmission_options.update({"WavelengthHigh": wav_high})
        if wav_step:
            transmission_options.update({"WavelengthStep": wav_step})
        if wav_step_type:
            transmission_options.update({"WavelengthStepType": wav_step_type})
        if rebin_mode:
            transmission_options.update({"RebinMode": rebin_mode})
        if fit_method:
            transmission_options.update({"FitMethod": fit_method})
        if poly_order:
            transmission_options.update({"PolynomialOrder": poly_order})
        if fit_wav_low:
            transmission_options.update({"FitWavelengthLow": fit_wav_low})
        if fit_wav_high:
            transmission_options.update({"FitWavelengthHigh": fit_wav_high})
        if fit_wav_step:
            transmission_options.update({"FitWavelengthStep": fit_wav_step})
        if fit_wav_step_type:
            transmission_options.update({"FitWavelengthStepType": fit_wav_step_type})
        transmisson_alg = create_unmanaged_algorithm(transmission_name, **transmission_options)
        transmisson_alg.execute()
        fitted = transmisson_alg.getProperty(SANSConstants.output_workspace).value
        unfitted = transmisson_alg.getProperty("Unfitted").value
        return fitted, unfitted

    def test_that(self):
        # Arrange
        transmission_workspace= None
        direct_workspace = None
        # Act
        self._perform_calculation()

        # Assert
        pass


class SANSCalculateTransmissionRunnerTest(stresstesting.MantidStressTest):
    def __init__(self):
        stresstesting.MantidStressTest.__init__(self)
        self._success = False

    def runTest(self):
        suite = unittest.TestSuite()
        suite.addTest(unittest.makeSuite(SANSCalculateTransmissionTest, 'test'))
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
