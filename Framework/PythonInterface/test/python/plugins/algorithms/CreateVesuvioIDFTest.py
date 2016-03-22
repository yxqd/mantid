import unittest
import math
from mantid.simpleapi import *
from IndirectImport import *
from mantid.api import MatrixWorkspace, WorkspaceGroup

class CreateVesuvioIDFTest(unittest.TestCase):

    _ip_filename = 'Vesuvio_IP_file_test.par'
    _idf_filename = 'generated_idf.xml'
    _old_definition = None
    _new_definition = None

    def tearDown(self):
        """
        Remove workspaces
        """
        if mtd.doesExist('old_idf'):
            mtd.remove('old_idf')
        if mtd.doesExist('new_idf'):
            mtd.remove('new_idf')

# --------------------------------- success cases ---------------------------------

    def test_generated_values(self):
        """
        Test LoadEmptyInstrument with generated IDF from IP file aganist LoadEmptyVesuvio with the same IP File
        Test will compare z values
        Due to current lack of Phi value, x and y can not be checked as they depend on Phi)
        """
        # Run new algorithm
        CreateVesuvioIDF(IPFilename=self._ip_filename,
                         OutputFilename=self._idf_filename)
        # Get Definitions
        LoadEmptyVesuvio(InstrumentParFile=self._ip_filename, OutputWorkspace='old_idf')
        LoadEmptyInstrument(Filename=self._idf_filename, OutputWorkspace='new_idf')
        self._old_definition = mtd['old_idf']
        self._new_definition = mtd['new_idf']

        # Obtain sample position and detector number for old/new files
        new_instrument = self._new_definition.getInstrument()
        old_instrument = self._old_definition.getInstrument()
        new_length = self._new_definition.getNumberHistograms()
        old_length = self._old_definition.getNumberHistograms()


        # Check number of detectors are equal
        self.assertEquals(new_length, old_length)

        # Loop over all detectors
        for i in range(0, new_length):
            # Get new/old detectors
            new_det = self._new_definition.getDetector(i)
            old_det = self._old_definition.getDetector(i)
            # Check for equality of z postions
            self.assertEquals(round(new_det.getPos().Z(), 5), round(old_det.getPos().Z(), 5))
            # Check for equality of IDs
            self.assertEquals(new_det.getID(), old_det.getID())
            # Check R values
            new_r = new_det.getDistance(new_instrument)
            old_r = old_det.getDistance(old_instrument)
            self.assertEquals(round(new_r,5), round(old_r,5))
            # Check Theta values
            new_theta = round((self._new_definition.detectorTwoTheta(new_det) * 180) / math.pi, 4)
            old_theta = round((self._old_definition.detectorTwoTheta(old_det) * 180) / math.pi, 4)
            self.assertEquals(new_theta, old_theta)


    def test_monitors(self):
        """
        Test that the monitors are 0 and 1
        """
        # Run new algorithm
        CreateVesuvioIDF(IPFilename=self._ip_filename,
                         OutputFilename=self._idf_filename)
        # Get Definitions
        LoadEmptyInstrument(Filename=self._idf_filename, OutputWorkspace='new_idf')
        self._new_definition = mtd['new_idf']

        # Get monitors
        monitor_one = self._new_definition.getDetector(0)
        monitor_two = self._new_definition.getDetector(1)

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
        self.assertEquals(monitor_two.getPhi(), 0)


#------------------------ failure cases -------------------------------------

    def test_invalid_file_extension(self):
        self.assertRaises(RuntimeError, CreateVesuvioIDF, self._ip_filename, 'wrong_extension.txt')


    def test_incorrect_file_length(self):
        ip_file = open('incorrect_length_vesuvio_ip_file.par', 'w')
        ip_file.write('single line')
        ip_file.close()

        self.assertRaises(RuntimeError, CreateVesuvioIDF, 'incorrect_length_vesuvio_ip_file.par', self._idf_filename)

        os.remove('incorrect_length_vesuvio_ip_file.par')


if __name__=="__main__":
    unittest.main()
