import unittest
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
        # Check instrument
        self.assertEquals(self._new_definition.getFullName(), 'VESUVIO')
        # Check number of detectors are equal
        new_length = self._new_definition.getNumberHistograms()
        old_length = self._old_definition.getNumberHistograms()
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
