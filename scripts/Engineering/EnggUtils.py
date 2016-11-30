#pylint: disable=invalid-name
#pylint: disable=too-many-arguments
# R0913: Too many arguments - strictly for write_ENGINX_GSAS_iparam_file()
from __future__ import (absolute_import, division, print_function)
from mantid.api import *
import mantid.simpleapi as sapi

INSTRUMENTS = ["ENGINX", "IMAT"]

ENGINX_BANKS = ['', 'North', 'South', 'Both: North, South', '1', '2']
ENGINX_MASK_BIN_MINS = [0, 19930, 39960, 59850, 79930]
ENGINX_MASK_BIN_MAXS = [5300, 20400, 40450, 62000, 82670]


def default_ceria_expected_peaks():
    """
    Get the list of expected Ceria peaks, which can be a good default for the expected peaks
    properties of algorithms like EnggCalibrate and EnggCalibrateFull

    @Returns :: a list of peaks in d-spacing as a float list
    """
    _CERIA_EXPECTED_PEAKS = [
        3.124277511, 2.705702376, 1.913220892, 1.631600313, 1.562138267,
        1.352851554, 1.241461538, 1.210027059, 1.104598643, 1.04142562,
        0.956610446, 0.914694494, 0.901900955, 0.855618487, 0.825231622,
        0.815800156, 0.781069134, 0.757748432, 0.750426918, 0.723129589,
        0.704504971, 0.676425777, 0.66110842, 0.656229382, 0.637740216,
        0.624855346, 0.620730846, 0.605013529
    ]

    return _CERIA_EXPECTED_PEAKS


def read_in_expected_peaks(filename, expectedGiven):
    """
    Reads in expected peaks from the .csv file if requested. Otherwise fall back to the list of
    peaks given (and check that it is not empty).

    @param :: filename name of the csv file to read from. If empty, we take the peaks given in an option.
    This is passed to Engg algorithms in the (optional) input property 'ExpectedPeaksFromFile'

    @param :: expectedGiven list of expected peaks given as an input property to an Engg algorithm
    (ExpectedPeaks)

    @returns the expected peaks either from a file or input list, sorted in ascending order
    """

    expectedPeaksD = None

    if filename:
        exPeakArray = []
        readInArray = []
        try:
            with open(filename) as f:
                for line in f:
                    readInArray.append([float(x) for x in line.split(',')])
            for a in readInArray:
                for b in a:
                    exPeakArray.append(b)
        except RuntimeError as exc:
            raise RuntimeError(
                "Error while reading file of expected peaks '%s': %s" %
                (filename, exc))

        if not exPeakArray:
            # "File could not be read. Defaults in alternative option used."
            if not expectedGiven:
                raise ValueError(
                    "Could not read any peaks from the file given in 'ExpectedPeaksFromFile: '"
                    + filename +
                    "', and no expected peaks were given in the property "
                    "'ExpectedPeaks' either. Cannot continue without a list of expected peaks."
                )

            expectedPeaksD = sorted(expectedGiven)

        else:
            expectedPeaksD = sorted(exPeakArray)

    else:
        if 0 == len(expectedGiven):
            raise ValueError(
                "No expected peaks were given in the property 'ExpectedPeaks', "
                "could not get default expected peaks, and 'ExpectedPeaksFromFile' "
                "was not given either. Cannot continout without a list of expected peaks."
            )

        expectedPeaksD = sorted(expectedGiven)

    return expectedPeaksD


def getWsIndicesFromInProperties(ws, bank, detIndices, switch):
    """
    Get the detector indices that the user requests, either through the input property 'Bank' or
    'DetectorIndices'

    @param workspace :: input workspace (with instrument)
    @param bank :: value passed in the input property 'Bank' to an Engg algorithm
    @param detIndices :: value passed in the 'Det

    @returns list of workspace indices that can be used in mantid algorithms such as CropWorkspace.
    """

    indices = None
    if bank and detIndices:
        raise ValueError(
            "It is not possible to use at the same time the input properties 'Bank' and "
            "'DetectorIndices', as they overlap. Please use either of them. Got Bank: '%s', "
            "and DetectorIndices: '%s'" % (bank, detIndices))

    elif bank:
        bankAliases = {'North': '1', 'South': '2', 'Both: North, South': '-1'}
        bank = bankAliases.get(bank, bank)
        indices = getWsIndicesForBank(ws, bank, switch)
        if not indices:
            raise RuntimeError(
                "Unable to find a meaningful list of workspace indices for the "
                "bank passed: %s. Please check the inputs." % bank)

        return indices
    elif detIndices:
        indices = parseSpectrumIndices(ws, detIndices)
        if not indices:
            raise RuntimeError(
                "Unable to find a meaningful list of workspace indices for the "
                "range(s) of detectors passed: %s. Please check the inputs." %
                detIndices)

        return indices
    else:
        raise ValueError(
            "You have not given any value for the properties 'Bank' and 'DetectorIndices' "
            "One of them is required")


def parseSpectrumIndices(ws, specNumbers):
    """
    Get a usable list of workspace indices from a user provided list of spectra that can look like:
    '8-10, 20-40, 100-110'. For that example this method will produce: [7,8,9,19, 20,... , 109]

    @param workspace :: input workspace (with instrument)
    @param specNumbers :: range of spectrum numbers (or list of ranges) as given to an algorithm

    @return list of workspace indices, ready to be used in mantid algorithms such as CropWorkspace
    """
    segments = [s.split("-") for s in specNumbers.split(",")]

    indices = [
        idx for s in segments for idx in range(int(s[0]), int(s[-1]) + 1)
    ]

    # remove duplicates and sort
    indices = sorted(set(indices))
    maxIdx = ws.getNumberHistograms()
    if indices[-1] >= maxIdx:
        raise ValueError(
            "A workspace index equal or bigger than the number of histograms available in the "
            "workspace '" + ws.getName() + "' (" + str(ws.getNumberHistograms(
            )) + ") has been given. Please check the list of indices.")

    # and finally traslate from 'spectrum numbers' to 'workspace indices'
    return [ws.getIndexFromSpectrumNumber(sn) for sn in indices]


def getBanksAndWSIndicesFor(ws, instrument=None):
    """
    Find the number of banks for the specified instrument, 
    and the corresponding workspace indices for each bank, 
    using the instrument's grouping file.

    ws :: workspace with instrument definition
    instrument :: instrument for which the grouping file will be read

    @returns banksDetectorIds :: list of banks containing list of indices, 
                each containing the workspace indices for its bank,
                the len(banksDetectorIds) will give the number of available banks
    """

    # Private funciton to determine if the index is in the bank
    banksDetectorIds = getBanksAndDetectorIDsFor(instrument)

    def isIndexInBank(index, bank):
        try:
            # we get the detector for each spectra, and check if that 
            # detector is in the list of all detectors that 
            # we loaded into banksDetectorIds
            det = ws.getDetector(index)
            return det.getID() in banksDetectorIds[bank]
        except RuntimeError:
            return False

    for b in range(len(banksDetectorIds)):
        # convert the detector IDs to their corresponding WS indices
        banksDetectorIds[b] = [
            i for i in range(0, ws.getNumberHistograms())
            if isIndexInBank(i, b)
        ]

    # at this point this contains the banks and their workspace indices
    return banksDetectorIds


def getBanksAndDetectorIDsFor(instrument=None,
                              fileSuffixAppend="_Grouping.xml"):
    """
    Find the detector IDs for an instrument bank. Note this is at this point specific to
    the ENGINX instrument.

    @param instrument :: name as a string.
    @param fileSuffixAppend :: specify the file suffix that will be 
            appended to the instrument name. The default is "_Grouping.xml"

    @returns banksDetectorIds :: list of banks containing list of detector IDs, 
            each containing the DETECTOR IDs for its bank,
            the len(banksDetectorIds) will give the number of available banks
    """
    import os
    filename = instrument + fileSuffixAppend

    groupingFilePath = os.path.join(sapi.config.getInstrumentDirectory(),
                                    'Grouping', filename)

    alg = AlgorithmManager.create('LoadDetectorsGroupingFile')
    alg.initialize()
    alg.setLogging(False)
    alg.setProperty('InputFile', groupingFilePath)
    grpName = '__IMATGrouping'
    alg.setProperty('OutputWorkspace', grpName)
    alg.execute()

    # LoadDetectorsGroupingFile produces a 'Grouping' workspace.
    # PropertyWithValue<GroupingWorkspace> not working (GitHub issue 13437)
    # => cannot run as child and get outputworkspace property properly
    if not AnalysisDataService.doesExist(grpName):
        raise RuntimeError(
            'LoadDetectorsGroupingFile did not run correctly. Could not '
            'find its output workspace: ' + grpName)

    grouping = mtd[grpName]

    banks = getDetectorIdsFromGroupings(grouping)

    sapi.DeleteWorkspace(grouping)

    if len(banks) == 0:
        raise ValueError('Could not find any detector for this bank: ' + bank +
                         '. This looks like an unknown bank')

    return banks


def getDetectorIdsFromGroupings(grouping):
    """
    This assumes that the banks numbers are in increasing order, 
    i.e. 1,2,3,4, in the groupings file, if they are not this will not work.
    
    @param groupings :: The groupings workspace produced by 
            the Mantid load groupings algorithm
    
    @returns banksDetectorIds :: list of banks containing list of detector IDs, 
        each containing the DETECTOR IDs for its bank,
        the len(banksDetectorIds) will give the number of available banks
    """

    banks = []
    banklen = len(banks)
    for i in range(grouping.getNumberHistograms()):
        # cast to integer, we cannot have bank 1.5
        currentDetectorBank = int(grouping.readY(i)[0])

        #exclude all detectors in bank 0 (i.e. excluded)
        if currentDetectorBank > 0:
            # append another bank if found
            if currentDetectorBank > banklen:
                banks.append([])
                # update bank len, so we only append for new banks
                banklen = len(banks)
            # we subtract 1 to get the bank list indice
            # this appends the actual detector ID
            banks[currentDetectorBank - 1].append(
                grouping.getDetector(i).getID())
    return banks


def getDetectorIdsFromAnyOrderInGroupings(groupings):
    """
    Not used, solves the problem with the banks not being in increasing order in the groupings file
    
    @param groupings :: The groupings workspace produced by 
            the Mantid load groupings algorithm
    
    @returns banksDetectorIds :: list of banks containing list of detector IDs, 
        each containing the DETECTOR IDs for its bank,
        the len(banksDetectorIds) will give the number of available banks
    """

    banks = []
    banklen = len(banks)
    for i in range(grouping.getNumberHistograms()):
        # cast to integer, we cannot have bank 1.5
        currentDetectorBank = int(grouping.readY(i)[0])

        #exclude all detectors in bank 0 (i.e. excluded)
        if currentDetectorBank > 0:
            # append another bank if found
            if currentDetectorBank > banklen:
                # this will append any number of difference in banks between
                # the current banks that we have (banklen) and the currentDetectorBank
                # this means if the first bank is 3, it will append 3 banks,
                # if the next banks is 1 it will not append, but just add to bank 1
                for i in range(currentDetectorBank - banklen):
                    banks.append([])
                # update bank len, so we only append for new banks
                banklen = len(banks)
            # we subtract 1 to get the bank list indice
            # this appends the actual detector ID
            banks[currentDetectorBank - 1].append(
                grouping.getDetector(i).getID())
    return banks


def generateOutputParTable(name, difa, difc, tzero):
    """
    Produces a table workspace with the two fitted calibration parameters

    @param name :: the name to use for the table workspace that is created here
    @param difa :: DIFA calibration parameter (GSAS parameter)
    @param difc :: DIFC calibration parameter
    @param tzero :: TZERO calibration parameter
    """
    tbl = sapi.CreateEmptyTableWorkspace(OutputWorkspace=name)
    tbl.addColumn('double', 'DIFA')
    tbl.addColumn('double', 'DIFZ')
    tbl.addColumn('double', 'TZERO')
    tbl.addRow([float(difa), float(difc), float(tzero)])


def applyVanadiumCorrections(parent, ws, vanWS, vanIntegWS, vanCurvesWS,
                             instrument):
    """
    Apply the EnggVanadiumCorrections algorithm on the workspace given, by using the algorithm
    EnggVanadiumCorrections

    @param parent :: parent (Mantid) algorithm that wants to run this
    @param ws :: workspace to correct (modified in place)
    @param vanWS :: workspace with data from a Vanadium run
    @param vanIntegWS :: alternatively to vanWS, pre-calculated integration from Vanadium data
    @param vanIntegWS :: alternatively to vanWS, pre-calculated bank curves from Vanadium data
    @param instrument :: the instrument for which the vanadium corrections will be calculated
    """
    # if vanWS and vanWS.getNumberHistograms() < len(indices):
    #     raise ValueError(
    #         "Inconsistency in inputs: the Vanadium workspace has less spectra (%d) than "
    #         "the number of workspace indices to process (%d)" %
    #         (vanWS.getNumberHistograms(), len(indices)))

    # elif vanIntegWS and vanCurvesWS:
    #     # filter only indices from vanIntegWS (crop the table)
    #     tbl = sapi.CreateEmptyTableWorkspace(
    #         OutputWorkspace="__vanadium_integration_ws")

    #     tbl.addColumn('double', 'Spectra Integration')
    #     for i in indices:
    #         tbl.addRow([vanIntegWS.cell(i, 0)])

    #     vanIntegWS = tbl

    # These corrections rely on ToF<->Dspacing conversions, so they're done after the calibration step
    # sapi.EnggVanadiumCorrections(Workspace=ws, VanadiumWorkspace=vanWS,
    #                              IntegrationWorkspace=vanIntegWS,
    #                              CurvesWorkspace=vanCurvesWS)
    alg = parent.createChildAlgorithm('EnggVanadiumCorrections')
    if ws:
        alg.setProperty('Workspace', ws)
    if vanWS:
        alg.setProperty('VanadiumWorkspace', vanWS)
    if vanIntegWS:
        alg.setProperty('IntegrationWorkspace', vanIntegWS)
    if vanCurvesWS:
        alg.setProperty('CurvesWorkspace', vanCurvesWS)
    if instrument:
        alg.setProperty('Instrument', instrument)
    alg.execute()


def convertToDSpacing(parent, ws):
    """
    Converts a workspace to dSpacing using 'ConvertUnits' as a child algorithm.

    @param parent :: parent (Mantid) algorithm that wants to run this
    @param ws :: workspace (normally in ToF units) to convert (not modified)

    @returns workspace converted to d-spacing units
    """
    # A check to catch possible errors in an understandable way
    expectedDim = 'Time-of-flight'
    dimType = ws.getXDimension().getName()
    if expectedDim != dimType:
        raise ValueError(
            "This function expects a workspace with %s X dimension, but "
            "the X dimension of the input workspace is: '%s'. This is an internal logic "
            "error. " % (expectedDim, dimType))

    alg = parent.createChildAlgorithm('ConvertUnits')
    alg.setProperty('InputWorkspace', ws)
    alg.setProperty('Target', 'dSpacing')
    alg.setProperty('AlignBins', True)
    alg.execute()
    return alg.getProperty('OutputWorkspace').value


def convertToToF(parent, ws):
    """
    Converts workspace to Time-of-Flight using 'ConvertUnits' as a child algorithm.

    @param parent :: parent (Mantid) algorithm that wants to run this
    @param ws :: workspace (normally in d-spacing units) to convert to ToF

    @returns workspace with data converted to ToF units
    """
    alg = parent.createChildAlgorithm('ConvertUnits')
    alg.setProperty('InputWorkspace', ws)
    alg.setProperty('Target', 'TOF')
    alg.execute()
    return alg.getProperty('OutputWorkspace').value


def cropData(parent, ws, indices):
    """
    Produces a cropped workspace from the input workspace so that only
    data for the specified bank (given as a list of indices) is left.

    NB: This assumes spectra for a bank are consequent.

    @param parent :: parent (Mantid) algorithm that wants to run this
    @param ws :: workspace to crop (not modified in-place)
    @param indices :: workspace indices to keep in the workpace returned

    @returns cropped workspace, with only the spectra corresponding to the indices requested
    """
    # Leave only spectra between min and max
    alg = parent.createChildAlgorithm('CropWorkspace')
    alg.setProperty('InputWorkspace', ws)
    alg.setProperty('StartWorkspaceIndex', min(indices))
    alg.setProperty('EndWorkspaceIndex', max(indices))
    alg.execute()

    return alg.getProperty('OutputWorkspace').value


def sumSpectra(parent, ws):
    """
    Focuses/sums up all the spectra into a single one (calls the SumSpectra algorithm)

    @param parent :: parent (Mantid) algorithm that wants to run this
    @param ws :: workspace to sum up

    @return single-spectrum workspace resulting from the sum
    """
    alg = parent.createChildAlgorithm('SumSpectra')
    alg.setProperty('InputWorkspace', ws)
    alg.execute()

    return alg.getProperty('OutputWorkspace').value


def write_ENGINX_GSAS_iparam_file(output_file,
                                  difc,
                                  tzero,
                                  bank_names=None,
                                  ceria_run=241391,
                                  vanadium_run=236516,
                                  template_file=None):
    """
    Produces and writes an ENGIN-X instrument parameter file for GSAS
    (in the GSAS iparam format, as partially described in the GSAS
    manual). It simply uses a template (found in template_path) where
    some values are replaced with the values (difc, tzero) passed to
    this function. DIFA is fixed to 0.

    Possible extensions for the file are .par (used here as default),
    .prm, .ipar, etc.

    @param output_file :: name of the file where to write the output
    @param difc :: list of DIFC values, one per bank, to pass on to GSAS
                   (as produced by EnggCalibrate)
    @param tzero :: list of TZERO values, one per bank, to pass on to GSAS
                    (also from EnggCalibrate)
    @param ceria_run :: number of the ceria (CeO2) run used for this calibration.
                        this number goes in the file and should also be used to
                        name the file
    @param vanadium_run :: number of the vanadium (VNb) run used for this
                           calibration. This number goes in the file and should
                           also be used to name the file.
    @param template_file :: file to use as template (with relative or full path)

    @returns

    """
    if not isinstance(difc, list) or not isinstance(tzero, list):
        raise ValueError(
            "The parameters difc and tzero must be lists, with as many elements as "
            "banks")

    if len(difc) != len(tzero):
        raise ValueError(
            "The lengths of the difc and tzero lists must be the same")

    # Defaults for a "both banks" file
    if not template_file:
        template_file = 'template_ENGINX_241391_236516_North_and_South_banks.prm'
    import os
    template_file = os.path.join(os.path.dirname(__file__), template_file)
    if not bank_names:
        bank_names = ["North", "South"]

    temp_lines = []
    with open(template_file) as tf:
        temp_lines = tf.readlines()

    def replace_patterns(line, patterns, replacements):
        """
        If line starts with any of the strings passed in the list 'pattern', return the
        corresponding 'replacement'
        """
        for idx, pat in enumerate(patterns):
            if line[0:len(pat)] == pat:
                return replacements[idx]

        return line

    # need to replace two types of lines/patterns:
    # - instrument constants/parameters (ICONS)
    # - instrument calibration comment with run numbers (CALIB)
    output_lines = []
    for b_idx, _bank_name in enumerate(bank_names):
        patterns = [
            "INS  %d ICONS" %
            (b_idx + 1),  # bank calibration parameters: DIFC, DIFA, TZERO
            "INS    CALIB",  # calibration run numbers (Vanadium and Ceria)
            "INS    INCBM"  # A his file for open genie (with ceria run number in the name)
        ]

        difa = 0.0
        # the ljust(80) ensures a length of 80 characters for the lines (GSAS rules...)
        replacements = [
            ("INS  {0} ICONS  {1:.2f}    {2:.2f}    {3:.2f}".format(
                b_idx + 1, difc[b_idx], difa, tzero[b_idx])).ljust(80) + '\n',
            ("INS    CALIB   {0}   {1} ceo2".format(ceria_run, vanadium_run)
             ).ljust(80) + '\n',
            ("INS    INCBM  ob+mon_{0}_North_and_South_banks.his".format(
                ceria_run)).ljust(80) + '\n'
        ]

        output_lines = [
            replace_patterns(line, patterns, replacements)
            for line in temp_lines
        ]

    with open(output_file, 'w') as of:
        of.writelines(output_lines)
