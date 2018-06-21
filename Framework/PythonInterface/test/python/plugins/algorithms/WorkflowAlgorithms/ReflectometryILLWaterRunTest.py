# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

from mantid.api import mtd
import numpy.testing
from testhelpers import assertRaisesNothing, create_algorithm
import unittest


class ReflectometryILLWaterRunTest(unittest.TestCase):
    def testDefaultRunExecutesSuccessfully(self):
        outWSName = 'outWS'
        args = {
            'Run': 'ILL/D17/317370.nxs',
            'OutputWorkspace': outWSName,
            'rethrow': True,
            'child': True
        }
        alg = create_algorithm('ReflectometryILLWaterRun', **args)
        assertRaisesNothing(self, alg.execute)
        outWS = alg.getProperty('OutputWorkspace').value
        self.assertEquals(outWS.getAxis(0).getUnit().caption(), 'Time-of-flight')
        mtd.clear()

    def testCleanupOFF(self):
        # test if intermediate workspaces exist:
        # not tested: position tables
        # normalise_to_slits, normalise_to_monitor, '_normalised_to_time_','transposed_flat_background'
        workspace_name_suffix = ['_water_merged_files_', '_water_detector_workspace_', '_water_scaled_workspace_']
        outWSName = 'outWS'
        args = {
            'Run': 'ILL/D17/317370.nxs',
            'OutputWorkspace': outWSName,
            'ScaleFactor': 4.,
            'Cleanup': 'Cleanup OFF',
            'SubalgorithmLogging': 'Logging ON',
            'rethrow': True,
            'child': True
        }
        alg = create_algorithm('ReflectometryILLWaterRun', **args)
        assertRaisesNothing(self, alg.execute)
        for i in range(len(workspace_name_suffix)):
            wsName = outWSName + workspace_name_suffix[i]
            self.assertTrue(mtd.doesExist(wsName), wsName)
        mtd.clear()

    if __name__ == "__main__":
        unittest.main()
