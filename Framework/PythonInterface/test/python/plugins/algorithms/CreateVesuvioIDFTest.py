import unittest
import math
from mantid.simpleapi import *
from IndirectImport import *
from mantid.api import MatrixWorkspace, WorkspaceGroup

class CreateVesuvioIDFTest(unittest.TestCase):

    _ip_filename = 'path_to_ip.par'
    _idf_filename = 'path_to_idf.xml'
    _old_definition = None
    _new_definition = None

    def setUp(self):
        """
        Generate IDF from IP file.
        Load instruments from IP and IDF
        """
        # Run new algorithm
        CreateVesuvioIDF(IPFilename=self._ip_filename,
                         OutputFilename=self._idf_filename)
        # Get Definitions
        self._old_definition = LoadEmptyVesuvio(InstrumentParFile=self._ip_filename)
        self._new_definition = LoadEmptyInstrument(Filename=self._idf_filename)

    def test_generated_values(self):
        """
        Test LoadEmptyInstrument with generated IDF from IP file aganist LoadEmptyVesuvio with the same IP File
        Test will compare z values
        Due to current lack of Phi value, x and y can not be checked as they depend on Phi)
        """
        # Obtain sample position and detector number for old/new files
        new_sample_pos = self._new_definition.getInstrument()
        old_sample_pos = self._old_definition.getInstrument()
        new_length = self._new_definition.getNumberHistograms()
        old_length = self._old_definition.getNumberHistograms()

        # Check Instrument is correct
        self.assertEquals(self._new_definition.getFullName(), 'VESUVIO')
        # Check number of detectors are equal
        self.assertEquals(new_length, old_length)

        # Loop over all detectors
        for i in range(0, new_length):
            # Get new/old detectors
            new_det = self._new_definition.getDetector(i)
            old_det = self._old_definition.getDetector(i)
            # Check for equality of z postions
            self.assertEquals(round(new_det.getPos()[2], 5), round(old_det.getPos()[2], 5))
            # Check for equality of IDs
            self.assertEquals(new_det.getID(), old_det.getID())
            # Check R values
            new_r = new_det.getDistance(new_sample_pos)
            old_r = old_det.getDistance(old_sample_pos)
            self.assertEquals(new_r, old_r)
            # Check Theta values
            new_theta = ((self._new_definition.detectorTwoTheta(new_det) * 180) / math.pi)
            old_theta = ((self._old_definition.detectorTwoTheta(old_det) * 180) / math.pi)
            self.assertEquals(new_theta, old_theta)


    def test_monitors(self):
        """
        Test that the monitors are 0 and 1
        """
        # Get monitors
        monitor_one = self._empty_instrument.getDetector(0)
        monitor_two = self._empty_instrument.getDetector(1)

        # validator monitors
        # Check monitor
        self.assertTrue(monitor_one.isMonitor())
        self.assertTrue(monitor_two.isMonitor())
        # Check Name
        self.assertEquals('Monitor 1', monitor_one.getName())
        self.assertEquals('Monitor 2', monitor_two.getName())
        # Check Phi
        self.assertEquals(monitor_one.getPhi(), 0)
        self.assertEquals(monitor_two.getPhi(), 0)


    if __name__=="__main__":
        unittest.main()
