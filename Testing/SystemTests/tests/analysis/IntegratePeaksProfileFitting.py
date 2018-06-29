#pylint: disable=no-init
from __future__ import (absolute_import, division, print_function)
import stresstesting
from mantid.simpleapi import *


class IntegratePeaksProfileFittingTest(stresstesting.MantidStressTest):
    '''Tests the IntegratePeaksProfileFitting workflow algorithm'''

    def runTest(self):
        LoadMD(Filename='MANDI_5921_ProfileFitting.nxs', LoadHistory=False, OutputWorkspace='profile_fitting_test_md')
        LoadIsawPeaks(Filename='MANDI_5921_ProfileFitting.integrate', OutputWorkspace='peaks_ws')
        IntegratePeaksProfileFitting(OutputPeaksWorkspace='peaks_ws_out', OutputParamsWorkspace='params_ws', 
                                     InputWorkspace='profile_fitting_test_md', PeaksWorkspace='peaks_ws', 
                                     RunNumber=5921, UBFile='/SNS/MANDI/shared/ProfileFitting/demo_5921.mat',
                                     ModeratorCoefficientsFile='/SNS/MANDI/shared/ProfileFitting/franz_coefficients_2017.dat', 
                                     StrongPeakParamsFile='/SNS/MANDI/shared/ProfileFitting/strongPeakParams_beta_lac_mut_mbvg.pkl', 
                                     MinpplFrac=0.99, MaxpplFrac=1.01, DtSpread = 0.015)

        intens3d = mtd['params_ws'].row(0)['Intens3d']
        sigint3d = mtd['params_ws'].row(0)['SigInt3d']
        self.assertDelta(intens3d, 281.07, 5., "Incorrect intensity found.")
        self.assertDelta(sigint3d, 20.168399, 1., "Incorrect sigma found.")

    def validate(self):
        return True
