# pylint: disable=no-init,invalid-name,attribute-defined-outside-init
from mantid.simpleapi import *
from mantid.api import *
from mantid.kernel import *

import numpy as np
from scipy.optimize import brent, fmin


def optimizationWrapperT0(t0, parameters, workspaces, algorithmObject):
    if np.fabs(t0) > 0.1:
        return 1e10

    paramCopy = [x for x in parameters]
    paramCopy[0] = t0

    # Slope differences are relevant
    # try:
    slopes = algorithmObject.getDataWithTimingParameters(workspaces, paramCopy)

    slopeDifferences = []
    for i in range(len(slopes)):
        for j in range(i + 1, len(slopes)):
            slopeDifferences.append(slopes[i] - slopes[j])

    return np.sum(np.square(np.array(slopeDifferences)))
    # except Exception as e:
    #    return 1e10


def optimizationWrapperTConst(tconst, parameters, t0, workspaces, algorithmObject):
    paramCopy = [x for x in parameters]
    paramCopy[0] = t0
    paramCopy[1] = tconst

    # Absolute values of slopes are checked for this parameter
    try:
        slopes = algorithmObject.getDataWithTimingParameters(workspaces, paramCopy)

        return np.sum(np.square(np.array(slopes)))
    except:
        return 1e10


class PoldiCalibration(PythonAlgorithm):
    """
    This workflow algorithm uses all of the POLDI specific algorithms to perform a complete data analysis,
    starting from the correlation method and preliminary 1D-fits, proceeding with either one or two passses
    of 2D-fitting.

    All resulting workspaces are grouped together at the end so that they are all in one place.
    """

    def category(self):
        return "SINQ\\Poldi"

    def name(self):
        return "PoldiCalibration"

    def summary(self):
        return "Calibrate POLDI using 5 parameters."

    def checkGroups(self):
        return False

    def PyInit(self):
        self.declareProperty(WorkspaceProperty(name='InputWorkspace', defaultValue='', direction=Direction.Input),
                             doc='WorkspaceGroup with POLDI runs at different chopper speeds.')

        self.declareProperty(WorkspaceProperty(name='ExpectedPeaks', defaultValue='', direction=Direction.Input),
                             doc='TableWorkspace with expected reflections')

        self.declareProperty('MaximumPeakNumber', 11, doc='Number of peaks to be used for calibration.')

        calibrationModes = StringListValidator(['Timing', 'Position'])
        self.declareProperty('CalibrationMode', 'Timing', validator=calibrationModes,
                             doc='Select which parameters are calibrated.')

        self.declareProperty('InitialParameters', '', direction=Direction.Input,
                             doc='Initial parameters for the calibration, comma separated in the order t0, '
                                 'tconst, x0, y0, two_theta. If not supplied, the loaded instrument parameters are '
                                 'used. For position refinement it is enough to provide the first two.')

        self.declareProperty('ParameterRanges', '', direction=Direction.Input,
                             doc='Parameter ranges for two_theta, x0, y0 in the form \'start,stop,number of points\', '
                                 'separated by semicolons. Ignored for t0, tconst optimization.')

        self.declareProperty(WorkspaceProperty(name='OutputWorkspace',
                                               defaultValue='', direction=Direction.Output),
                             doc='Output table with parameters for position calibration.')

    def PyExec(self):
        workspaceGroup = self.getProperty('InputWorkspace').value

        if not isinstance(workspaceGroup, WorkspaceGroup):
            raise RuntimeError('InputWorkspace needs to be a WorkspaceGroup.')

        # Get input workspaces
        workspaceList = [AnalysisDataService.retrieve(x) for x in workspaceGroup.getNames()]
        calibrationMode = self.getProperty('CalibrationMode').value

        # Get expected peaks
        self._expectedPeaks = self.getProperty('ExpectedPeaks').value
        self._maxPeakCount = self.getProperty('MaximumPeakNumber').value

        # Store initial parameters
        initialParameters = self.getProperty('InitialParameters')
        if not initialParameters.isDefault:
            self.setInitialParametersFromString(initialParameters.value)
        else:
            self.setInitialParametersFromWorkspace(workspaceList[0])

        # Perform calibration
        if calibrationMode == 'Timing':
            t0, tconst = self.calibrateTiming(workspaceList)

            self.log().warning('''Calibrated timing parameters:
                t0 = {}
                tconst = {}'''.format(t0, tconst))

            outputWorkspace = self.createOutputWorkspaceTiming(t0, tconst)
            self.setProperty('OutputWorkspace', outputWorkspace)
        else:
            lines = self.calibratePosition(workspaceList)

            outputWorkspace = self.createOutputWorkspacePosition(lines)
            self.setProperty('OutputWorkspace', outputWorkspace)

    def createOutputWorkspacePosition(self, data):
        columnNames = ['TwoTheta', 'x0', 'y0', 'a', 'delta_a', 'slope', 'delta_slope', 'fwhm', 'delta_fwhm']

        tableWs = WorkspaceFactory.createTable()
        for n in columnNames:
            tableWs.addColumn('double', n)

        for row in data:
            tableWs.addRow(dict(zip(columnNames, row)))

        return tableWs

    def createOutputWorkspaceTiming(self, t0, tconst):
        tableWs = WorkspaceFactory.createTable()
        tableWs.addColumn('str', 'Parameter')
        tableWs.addColumn('double', 'Value')

        tableWs.addRow({'Parameter': 't0', 'Value': t0})
        tableWs.addRow({'Parameter': 'tconst', 'Value': tconst})

        return tableWs

    def setInitialParametersFromWorkspace(self, workspace):
        instrument = workspace.getInstrument()

        chopper = instrument.getComponentByName('chopper')
        t0 = chopper.getNumberParameter('t0')[0]
        tconst = chopper.getNumberParameter('t0_const')[0]

        detector = instrument.getComponentByName('detector')
        position = detector.getPos()
        two_theta = detector.getNumberParameter('two_theta')[0]

        self._initialParameters = [t0, tconst, position.X(), position.Y(), two_theta]

    def setInitialParametersFromString(self, parameterString):
        parts = [float(x) for x in parameterString.split(',')]

        if len(parts) == 2:
            parts += [0.0] * 3

        if not len(parts) == 5:
            raise RuntimeError('InitialParameters format is wrong.')

        self._initialParameters = parts

    def calibratePosition(self, workspaces):
        workspace = self.getWorkspaceWithHighestChopperSpeed(workspaces)
        ranges = self.getRangesFromProperty()

        combinations = np.cumprod([len(x) for x in ranges])[-1]
        self.log().warning('Number of parameter combinations: ' + str(combinations))

        lines = []

        prog = Progress(self, 0.0, 1.0, combinations)

        for two_theta in ranges[0]:
            for x0 in ranges[1]:
                for y0 in ranges[2]:
                    params = self.getParameters(two_theta, x0, y0)

                    self.log().warning('Parameters: ' + str(params))

                    ws = self.getWorkspaceWithParameters(workspace, *params)

                    slope, slope_error = self.getSlopeParameter(ws)
                    a, a_error, fwhms = self.getLatticeParameter(ws)

                    lines.append([two_theta, x0, y0, a, a_error, slope * 1000.0, slope_error * 1000.0,
                                  fwhms[0][0], fwhms[0][1]])

                    prog.report()

        return lines

    def saveDataPoints(self, datalines):
        # Expected a list of lists of floats
        fileName = self.getProperty('Output').value

        fh = open(fileName, 'w')
        fh.write('# two_theta x0 y0 a delta_a slope delta_slope fwhm delta_fwhm')
        for l in datalines:
            fh.write(' '.join(l))
            fh.write('\n')

        fh.close()

    def getParameters(self, two_theta, x0, y0):
        params = [x for x in self._initialParameters[:2]]
        params += [x0, y0, two_theta]

        return params

    def _getChopperSpeed(self, ws):
        try:
            return ws.getRun().getProperty('chopperspeed').value[0]
        except Exception:
            return ws.getRun().getProperty('chopperspeed').value

    def getWorkspaceWithHighestChopperSpeed(self, workspaces):
        resultWs = workspaces[0]

        highestChopperSpeed = self._getChopperSpeed(resultWs)

        for ws in workspaces[1:]:
            chopperSpeed = self._getChopperSpeed(ws)
            print chopperSpeed

            if chopperSpeed > highestChopperSpeed:
                highestChopperSpeed = chopperSpeed
                resultWs = ws

        return resultWs

    def getRangesFromProperty(self):
        rangeString = self.getProperty('ParameterRanges').value
        rangeStrings = rangeString.split(';')

        if not len(rangeStrings) == 3:
            raise RuntimeError('Three parameter ranges must be specified.')

        ranges = []
        for rStr in rangeStrings:
            params = [float(x) for x in rStr.split(',')]
            if not len(params) == 3:
                raise RuntimeError('Parameter range must be specified by start,stop,number of points')

            ranges.append(np.linspace(start=min(params[:2]), stop=max(params[:2]), num=int(params[2])))

        return ranges

    def calibrateTiming(self, workspaces):
        # First, calibrate t0
        t0 = self.calibrateT0(workspaces)

        # Since a rounded value will be put into the parameter file, this is done here, both values are logged
        self.log().notice('Calibrated value for t0: ' + str(t0))
        t0_rounded = np.round(t0, 6)
        self.log().notice('Rounded value for t0 used as parameter: ' + str(t0_rounded))

        # With the new t0 value, tconst can be refined as well.
        tconst = self.calibrateTConst(workspaces, t0_rounded)
        self.log().notice('Calibrated value for tconst: ' + str(tconst))
        tconst_rounded = np.round(tconst, 3)
        self.log().notice('Rounded value for tconst used as parameter: ' + str(tconst_rounded))

        return t0_rounded, tconst_rounded

    def calibrateT0(self, workspaces):
        return brent(optimizationWrapperT0, args=(self._initialParameters, workspaces, self), brack=(-0.09, 0.01),
                     tol=1e-4)

    def calibrateTConst(self, workspaces, t0):
        return brent(optimizationWrapperTConst, args=(self._initialParameters, t0, workspaces, self),
                     brack=(-20.0, 20.0), tol=1e-3)

    def getDataWithTimingParameters(self, workspaces, parameters):
        self.log().warning('Parameters: ' + str(parameters))

        slopes = []
        for ws in workspaces:
            realWs = self.getWorkspaceWithParameters(ws, *parameters)
            slopes.append(self.getLatticeParameterSlope(realWs))

        return slopes

    def getLatticeParameterSlope(self, workspace):
        fitResult = PoldiDataAnalysis(InputWorkspace=workspace, ExpectedPeaks=self._expectedPeaks,
                                      MaximumPeakNumber=self._maxPeakCount, AnalyseResiduals=False)

        # fitResult is a workspaceGroup and the refined peaks are in the last workspace.
        fittedPeaks = fitResult.getItem(fitResult.size() - 1)

        # if there are unindexed peaks, this is a workspace group
        if isinstance(fittedPeaks, WorkspaceGroup):
            fittedPeaks = fittedPeaks.getItem(0)

        # Extract hkls and q-values
        hkls = fittedPeaks.column(0)

        q_values = np.array(fittedPeaks.column(3))
        q_errors = np.array(fittedPeaks.column(4))
        rel_errors = q_errors / q_values

        # calculate sqrt(h^2 + k^2 + l^2) - calibration substance is always cubic
        hkl_sums = [np.sqrt(np.sum([int(x) * int(x) for x in y.split()])) for y in hkls]

        # calculate a values and errors
        a_values = np.array([y * (2.0 * np.pi / x) for x, y in zip(q_values, hkl_sums)])
        a_errors = np.array([x * y for x, y in zip(rel_errors, a_values)])

        # remove outliers
        median = np.median(a_values)
        iqr = np.percentile(a_values, 75) - np.percentile(a_values, 25)

        goodValues = np.fabs(a_values - median) < 4.0 * iqr

        a_values_good = a_values[goodValues]
        a_errors_good = a_errors[goodValues]
        q_values_good = q_values[goodValues]

        self.log().notice('Number of peaks used for analysis: ' + str(len(a_values_good)))

        # fit a linear function to the data points
        slopeWs = CreateWorkspace(q_values, a_values, a_errors)

        fitStatus, chiSq, covarianceTable, paramTable, fitWorkspace = Fit("name=LinearBackground",
                                                                          IgnoreInvalidData=True,
                                                                          InputWorkspace=slopeWs, CreateOutput=True)
        slope = paramTable.cell(1, 1)

        self.log().information('Slope: ' + str(slope))
        self.log().information('Chi^2 of fit: ' + str(chiSq))

        covarianceTable.delete()
        paramTable.delete()
        # fitWorkspace.delete()
        fitResult.delete()

        return slope

    def getSlopeParameter(self, workspace):
        # Do a calibration run, which computes the additional slope parameter
        # Note: Despite the similarity in name, this has nothing to do with the slope in lattice parameters above.
        fitResult = PoldiDataAnalysis(InputWorkspace=workspace, CalibrationRun=True, OutputRawFitParameters=True,
                                      ExpectedPeaks=self._expectedPeaks, MaximumPeakNumber=self._maxPeakCount,
                                      AnalyseResiduals=False)

        # Get the slope parameter and return it.
        rawParams = fitResult.getItem(fitResult.size() - 1)
        slope = float(rawParams.cell(3, 1))
        slope_error = float(rawParams.cell(3, 2))

        fitResult.delete()

        return slope, slope_error

    def getLatticeParameter(self, workspace):
        fitResult = PoldiDataAnalysis(InputWorkspace=workspace, OutputRawFitParameters=True,
                                      ExpectedPeaks=self._expectedPeaks,
                                      MaximumPeakNumber=self._maxPeakCount, AnalyseResiduals=False)
        fpeaks = AnalysisDataService.retrieve(workspace.getName() + '_peaks_refined_2d')

        if isinstance(fpeaks, WorkspaceGroup):
            fpeaks = fpeaks.getItem(0)

        fwhms = zip(fpeaks.column(7), fpeaks.column(8))

        Fit("name=LatticeFunction,LatticeSystem=Cubic,a=5.4", Ties="ZeroShift=0.0",
            CostFunction="Unweighted least squares", InputWorkspace=fpeaks, CreateOutput=True, Output='cell')

        params = AnalysisDataService.retrieve('cell_Parameters')
        values = params.column(1)
        errors = params.column(2)

        return float(values[0]), float(errors[0]), fwhms

    def getWorkspaceWithParameters(self, workspace, t0, tconst, x0, y0, two_theta):
        workWs = workspace.clone()

        MoveInstrumentComponent(workWs, ComponentName='detector', x=x0, y=y0, RelativePosition=False)
        SetInstrumentParameter(workWs, ParameterName="two_theta", ComponentName="detector", ParameterType="Number",
                               Value=str(two_theta))

        SetInstrumentParameter(workWs, ParameterName="t0", ComponentName="chopper", ParameterType="Number",
                               Value=str(t0))
        SetInstrumentParameter(workWs, ParameterName="t0_const", ComponentName="chopper", ParameterType="Number",
                               Value=str(tconst))

        return workWs


AlgorithmFactory.subscribe(PoldiCalibration)
