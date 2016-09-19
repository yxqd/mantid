import unittest
import stresstesting

import mantid
from SANS2.Common.SANSConstants import SANSConstants
from SANS2.Common.SANSFunctions import create_unmanaged_algorithm


class SANSCalculateTransmissionTest(unittest.TestCase):
    @staticmethod
    def _load_data(file_name):
        load_name = "LoadISISNexus"
        load_options = {SANSConstants.file_name: file_name,
                        SANSConstants.output_workspace: SANSConstants.dummy}
        load_alg = create_unmanaged_algorithm(load_name, **load_options)
        load_alg.execute()
        return load_alg.getProperty(SANSConstants.output_workspace).value

    def _perform_calculation(self, trans_workspace, direct_workspace, incid_mon=None, trans_mon=None, trans_radius=None,
                             trans_roi=None, trans_mask=None, prompt_start=None, prompt_stop=None, flat_start=None,
                             flat_stop=None, flat_mon=None, wav_low=None, wav_high=None, wav_step=None,
                             wav_step_type=None, rebin_mode=None, fit_method=None, poly_order=None, fit_wav_low=None,
                             fit_wav_high=None, fit_wav_step=None, fit_wav_step_type=None):
        transmission_name = "SANSCalculateTransmission"
        transmission_options = {"TransmissionWorkspace": trans_workspace,
                                "DirectWorkspace": direct_workspace,
                                SANSConstants.output_workspace: SANSConstants.dummy}
        if incid_mon is not None:
            transmission_options.update({"IncidentMonitorSpectrumNumber": incid_mon})
        # --------------
        if trans_mon is not None:
            transmission_options.update({"TransmissionMonitor": trans_mon})
        if trans_radius is not None:
            transmission_options.update({"TransmissionRadius": trans_radius})
        if trans_roi is not None:
            transmission_options.update({"TransmissionROIFiles": trans_roi})
        if trans_mask is not None:
            transmission_options.update({"TransmissionMaskFiles": trans_mask})
        # --------------
        if prompt_start is not None:
            transmission_options.update({"PromptPeakCorrectionStart": prompt_start})
        if prompt_stop is not None:
            transmission_options.update({"PromptPeakCorrectionStop": prompt_stop})
        # --------------
        if flat_start is not None:
            transmission_options.update({"FlatBackgroundCorrectionROIStart": flat_start})
        if flat_stop is not None:
            transmission_options.update({"PromptPeakCorrectionStop": flat_stop})
        if flat_mon is not None:
            transmission_options.update({"FlatBackgroundCorrectionMonitors": flat_mon})
        # --------------
        if wav_low is not None:
            transmission_options.update({"WavelengthLow": wav_low})
        if wav_high is not None:
            transmission_options.update({"WavelengthHigh": wav_high})
        if wav_step is not None:
            transmission_options.update({"WavelengthStep": wav_step})
        if wav_step_type is not None:
            transmission_options.update({"WavelengthStepType": wav_step_type})
        if rebin_mode is not None:
            transmission_options.update({"RebinMode": rebin_mode})
        if fit_method is not None:
            transmission_options.update({"FitMethod": fit_method})
        if poly_order is not None:
            transmission_options.update({"PolynomialOrder": poly_order})
        if fit_wav_low is not None:
            transmission_options.update({"FitWavelengthLow": fit_wav_low})
        if fit_wav_high is not None:
            transmission_options.update({"FitWavelengthHigh": fit_wav_high})
        if fit_wav_step is not None:
            transmission_options.update({"FitWavelengthStep": fit_wav_step})
        if fit_wav_step_type is not None:
            transmission_options.update({"FitWavelengthStepType": fit_wav_step_type})

        print "=============================="
        print transmission_options

        transmission_alg = create_unmanaged_algorithm(transmission_name, **transmission_options)
        transmission_alg.execute()
        fitted = transmission_alg.getProperty(SANSConstants.output_workspace).value
        unfitted = transmission_alg.getProperty("Unfitted").value
        return fitted, unfitted

    def test_that_standard_transmission_monitor_calculation_can_be_performed(self):
        # Arrange
        direct_workspace = self._load_data("SANS2D00022024")
        transmission_workspace = self._load_data("SANS2D00022041")
        # Act
        flat_monitors = {"1": [80000, 90000], "2": [80000, 90000]}
        fitted, unfitted = self._perform_calculation(transmission_workspace, direct_workspace, incid_mon=1,
                                                     trans_mon=2, flat_start=80000., flat_stop=90000.,
                                                     flat_mon=flat_monitors, wav_low=4, wav_high=14.5, wav_step=2.,
                                                     wav_step_type="LIN",  rebin_mode="Rebin", fit_method="Linear",
                                                     poly_order=None, fit_wav_low=2., fit_wav_high=14.5,
                                                     fit_wav_step=2., fit_wav_step_type="LIN")
        # Assert
        print "========================="
        print "((((((((((((((((((((("
        print type(fitted)
        print type(unfitted)


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
