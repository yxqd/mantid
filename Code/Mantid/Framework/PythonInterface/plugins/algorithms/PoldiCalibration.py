# pylint: disable=no-init,invalid-name,attribute-defined-outside-init
from mantid.simpleapi import *
from mantid.api import *
from mantid.kernel import *

import numpy as np
from scipy.optimize import brent, fmin


def optimizationWrapperT0(t0, parameters, workspaces, algorithmObject):
    if np.fabs(t0[0]) > 0.1:
        return np.array([1.e10])

    paramCopy = [x for x in parameters]
    paramCopy[0] = t0[0]

    # Slope differences are relevant
    try:
        slopes = algorithmObject.getDataWithTimingParameters(workspaces, paramCopy)

        slopeDifferences = []
        for i in range(len(slopes)):
            for j in range(i + 1, len(slopes)):
                slopeDifferences.append(slopes[i] - slopes[j])

        return np.array([np.sum(np.square(np.array(slopeDifferences)))])
    except:
        return np.array([1.e10])


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
                                 'used.')

        self.declareProperty('ParameterRanges', '', direction=Direction.Input,
                             doc='Parameter ranges for two_theta, x0, y0 in the form \'start,stop,number of steps\', '
                                 'separated by semicolons. Ignored for t0, tconst optimization.')

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
            self.calibrateTiming(workspaceList)
        else:
            self.calibratePosition(workspaceList)

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

        if not len(parts) == 5:
            raise RuntimeError('InitialParameters format is wrong.')

        self._initialParameters = parts

    def calibratePosition(self, workspaces):
        pass


    def calibrateTiming(self, workspaces):
        #t0 = self.calibrateT0(workspaces)
        t0 = -0.003375

        self.log().warning('Calibrated value for t0: ' + str(t0))
        t0_rounded = np.round(t0, 6)
        self.log().warning('Rounded value for t0 used as parameter: ' + str(t0_rounded))

        tconst = self.calibrateTConst(workspaces, t0_rounded)
        self.log().warning('Calibrated value for tconst: ' + str(tconst))
        tconst_rounded = np.round(tconst, 3)
        self.log().warning('Rounded value for tconst used as parameter: ' + str(tconst_rounded))

        print t0_rounded, tconst_rounded


    def calibrateT0(self, workspaces):
        #return brent(optimizationWrapperT0, args=(self._initialParameters, workspaces, self), brack=(-0.09, 0.01),
        #             tol=1e-4)
        return fmin(optimizationWrapperT0, x0=np.array([0.0]), args=(self._initialParameters, workspaces,
                                                                              self),
                             xtol=1e-4, ftol=1e-8)[0]

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

        qs = fittedPeaks.column(2)
        q_values = np.array([float(x.split()[0]) for x in qs])
        q_errors = [float(x.split()[-1]) for x in qs]
        rel_errors = [x / y for x, y in zip(q_errors, q_values)]

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
        slopeWs = CreateWorkspace(q_values_good, a_values_good, a_errors_good)

        fitStatus, chiSq, covarianceTable, paramTable, fitWorkspace = Fit("name=LinearBackground",
                                                                          IgnoreInvalidData=True,
                                                                          InputWorkspace=slopeWs, CreateOutput=True)
        slope = paramTable.cell(1, 1)

        self.log().information('Slope: ' + str(slope))
        self.log().information('Chi^2 of fit: ' + str(chiSq))

        return slope

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


AlgorithmFactory.subscribe(PoldiCalibration())