# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

from mantid.api import AlgorithmFactory, DataProcessorAlgorithm, FileAction, MatrixWorkspaceProperty,\
    MultipleFileProperty, PropertyMode
from mantid.kernel import Direction, Property, StringListValidator
from mantid.simpleapi import ExtractMonitors, LoadAndMerge, mtd, Scale
import ReflectometryILL_common as common


class Prop():
    CLEANUP = 'Cleanup'
    OUTPUT_WS = 'OutputWorkspace'
    RUN = 'Run'
    SCALE_FACTOR = 'ScaleFactor'
    SUBALG_LOGGING = 'SubalgorithmLogging'


class SubalgLogging():
    OFF = 'Logging OFF'
    ON = 'Logging ON'


class ReflectometryILLWaterRun(DataProcessorAlgorithm):

    def category(self):
        """Return the categories of the algrithm."""
        return 'ILL\\Reflectometry;Workflow\\Reflectometry'

    def name(self):
        """Return the name of the algorithm."""
        return 'ReflectometryILLWaterRun'

    def summary(self):
        """Return a summary of the algorithm."""
        return "Loads, merges and normalizes ILL water file reflectometry data."

    def seeAlso(self):
        """Return a list of related algorithm names."""
        return ['ReflectometryILLConvertToQ', 'ReflectometryILLPolarizationCor', 'ReflectometryILLPreprocess',
                'ReflectometryILLSumForeground']

    def version(self):
        """Return the version of the algorithm."""
        return 1

    def PyExec(self):
        """Execute the algorithm."""
        self._subalgLogging = self.getProperty(Prop.SUBALG_LOGGING).value == SubalgLogging.ON
        cleanupMode = self.getProperty(Prop.CLEANUP).value
        self._cleanup = common.WSCleanup(cleanupMode, self._subalgLogging)
        wsPrefix = self.getPropertyValue(Prop.OUTPUT_WS)
        self._names = common.WSNameSource(wsPrefix, cleanupMode)

        run = self.getPropertyValue(Prop.RUN)

        mergedWSName = self._names.withSuffix('water_merged_files')
        ws = LoadAndMerge(
             Filename=run,
             LoaderName="LoadILLReflectometry",
             LoaderOptions={'XUnit' : 'TimeOfFlight'},
             OutputWorkspace=mergedWSName,
             EnableLogging=self._subalgLogging
        )

        detectorWSName = self._names.withSuffix('water_detector_workspace')
        ExtractMonitors(
            InputWorkspace=ws,
            DetectorWorkspace=detectorWSName,
            EnableLogging=self._subalgLogging
        )
        ws = mtd[detectorWSName]

        if not self.getProperty(Prop.SCALE_FACTOR).isDefault:
            scaledWSName = self._names.withSuffix('water_scaled_workspace')
            ws = Scale(
                 InputWorkspace=ws,
                 OutputWorkspace=scaledWSName,
                 Factor=self.getProperty(Prop.SCALE_FACTOR).value,
                 EnableLogging=self._subalgLogging
            )

        self._finalize(ws)

    def PyInit(self):
        """Initialize the input and output properties of the algorithm."""
        self.declareProperty(MultipleFileProperty(Prop.RUN,
                                                  action=FileAction.Load,
                                                  extensions=['nxs']),
                             doc='A list of water input run numbers/files.')
        self.declareProperty(MatrixWorkspaceProperty(Prop.OUTPUT_WS,
                                                     defaultValue='',
                                                     direction=Direction.Output),
                             doc='The preprocessed water run (unit TOF).')
        self.declareProperty(Prop.SCALE_FACTOR,
                             defaultValue=Property.EMPTY_DBL,
                             doc='Scale factor.')
        self.declareProperty(Prop.SUBALG_LOGGING,
                             defaultValue=SubalgLogging.OFF,
                             validator=StringListValidator([SubalgLogging.OFF, SubalgLogging.ON]),
                             doc='Enable or disable child algorithm logging.')
        self.declareProperty(Prop.CLEANUP,
                             defaultValue=common.WSCleanup.ON,
                             validator=StringListValidator([common.WSCleanup.ON, common.WSCleanup.OFF]),
                             doc='Enable or disable intermediate workspace cleanup.')

    def _finalize(self, ws):
        """Set OutputWorkspace to ws and clean up."""
        self.setProperty(Prop.OUTPUT_WS, ws)
        self._cleanup.cleanup(ws)
        self._cleanup.finalCleanup()


AlgorithmFactory.subscribe(ReflectometryILLWaterRun)
