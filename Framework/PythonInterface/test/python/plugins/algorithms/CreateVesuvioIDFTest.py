import unittest
from mantid.simpleapi import *
from IndirectImport import *
from mantid.api import MatrixWorkspace, WorkspaceGroup

class CreateVesuvioIDFTest(unittest.TestCase):

    _ip_filename = 'path_to_ip.par'
    _idf_filename = 'path_to_idf.xml'
    _empty_vesuvio = None
    _empty_instrument = None

    def setUp(self):
        """
        Generate IDF from IP file.
        Load instruments from IP and IDF
        """
        CreateVesuvioIDF(IPFilename=self._ip_filename,
                         OutputFilename=self._idf_filename)
        self._empty_vesuvio = LoadEmptyVesuvio(InstrumentParFile=self._ip_filename)
        self._empty_inst = LoadEmptyInstrument(Filename=self._idf_filename)

    def test_generated_idf(self):
        """
        Test LoadEmptyInstrument with generated IDF from IP file aganist LoadEmptyVesuvio with the same IP File
        """






    if __name__=="__main__":
        unittest.main()
