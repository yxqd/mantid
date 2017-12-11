from __future__ import (absolute_import, division, print_function)

from mantid.simpleapi import (AppendSpectra, CloneWorkspace, ElasticWindow, LoadLog, Logarithm, SortXAxis, Transpose)
from mantid.kernel import *
from mantid.api import *

import numpy as np


def _normalize_by_index(workspace, index):
    """
    Normalize each spectra of the specified workspace by the y-value at the
    specified index in that spectra, while accounting for errors on those values.

    @param workspace    The workspace to normalize.
    @param index        The index of the y-value to normalize by.
    """

    num_hist = workspace.getNumberHistograms()

    # Normalize each spectrum in the workspace
    for idx in range(0, num_hist):
        y_vals = workspace.readY(idx)
        y_vals_e = workspace.readE(idx)
        scale = y_vals[index]
        scale_e = y_vals_e[index]
        y_vals_scaled = y_vals / scale
        # error propagation: C = A / B ; dC = sqrt( (dA/B)^2 + (A*dB/B^2)^2 ) ||
        # !! wrong for A=B (index by which is scaled = index) !!
        if idx == index:
            y_vals_e_scaled = y_vals_e / scale
        else:
            y_vals_e_scaled = np.sqrt((y_vals_e / scale) ** 2 + (y_vals_scaled * scale_e / scale) ** 2)
        workspace.setY(idx, y_vals_scaled)
        workspace.setE(idx, y_vals_e_scaled)


def _append_workspaces(workspace1, workspace2):
    """
    Appends the spectra of the second workspace to the first.

    :param workspace1: The workspace to append spectra to.
    :param workspace2: The workspace whose spectra to append.
    :return:           The workspace containing the second workspace appended to the first.
    """
    return AppendSpectra(InputWorkspace1=workspace1, InputWorkspace2=workspace2,
                         OutputWorkspace="appended", StoreInADS=False, EnableLogging=False)

def _transpose_workspace(workspace):
    """
    Transposes the specified workspace.

    :param workspace:   The workspace to transpose.
    :return:            The tranposed workspace.
    """
    return Transpose(InputWorkspace=workspace, OutputWorkspace="transpose",
                     StoreInADS=False, EnableLogging=False)

def _clone_workspace(workspace):
    """
    Clones the specified workspace.

    :param workspace:   The workspace to clone.
    :return:            The cloned workspace.
    """
    return CloneWorkspace(InputWorkspace=workspace, OutputWorkspace="cloned",
                          StoreInADS=False, EnableLogging=False)

def _sort_workspace(workspace):
    """
    Sorts the specified workspace by it's X-Axis.

    :param workspace:   The workspace to sort.
    :return:            The sorted workspace.
    """
    return SortXAxis(InputWorkspace=workspace, OutputWorkspace="sorted",
                     StoreInADS=False, EnableLogging=False)


class ElasticWindowMultiple(DataProcessorAlgorithm):
    _sample_log_name = None
    _sample_log_value = None
    _input_workspaces = None
    _integration_range_start = None
    _integration_range_end = None
    _background_range_start = None
    _background_range_end = None

    def category(self):
        return 'Workflow\\Inelastic;Inelastic\\Indirect'

    def summary(self):
        return 'Performs the ElasticWindow algorithm over multiple input workspaces'

    def PyInit(self):
        self.declareProperty(WorkspaceGroupProperty('InputWorkspaces', '', Direction.Input),
                             doc='Grouped input workspaces')

        self.declareProperty(name='IntegrationRangeStart', defaultValue=0.0,
                             doc='Start of integration range in time of flight')
        self.declareProperty(name='IntegrationRangeEnd', defaultValue=0.0,
                             doc='End of integration range in time of flight')

        self.declareProperty(name='BackgroundRangeStart', defaultValue=Property.EMPTY_DBL,
                             doc='Start of background range in time of flight')
        self.declareProperty(name='BackgroundRangeEnd', defaultValue=Property.EMPTY_DBL,
                             doc='End of background range in time of flight')

        self.declareProperty(name='SampleEnvironmentLogName', defaultValue='sample',
                             doc='Name of the sample environment log entry')

        sampEnvLogVal_type = ['last_value', 'average']
        self.declareProperty('SampleEnvironmentLogValue', 'last_value',
                             StringListValidator(sampEnvLogVal_type),
                             doc='Value selection of the sample environment log entry')

        self.declareProperty(WorkspaceProperty('OutputInQ', '', Direction.Output),
                             doc='Output workspace in Q')

        self.declareProperty(WorkspaceProperty('OutputInQSquared', '', Direction.Output),
                             doc='Output workspace in Q Squared')

        self.declareProperty(WorkspaceProperty('OutputELF', '', Direction.Output,
                                               PropertyMode.Optional),
                             doc='Output workspace ELF')

        self.declareProperty(WorkspaceProperty('OutputELT', '', Direction.Output,
                                               PropertyMode.Optional),
                             doc='Output workspace ELT')

    def validateInputs(self):
        issues = dict()

        background_range_start = self.getProperty('BackgroundRangeStart').value
        background_range_end = self.getProperty('BackgroundRangeEnd').value

        if background_range_start != Property.EMPTY_DBL and background_range_end == Property.EMPTY_DBL:
            issues['BackgroundRangeEnd'] = 'If background range start was given and ' \
                                           'background range end must also be provided.'

        if background_range_start == Property.EMPTY_DBL and background_range_end != Property.EMPTY_DBL:
            issues['BackgroundRangeStart'] = 'If background range end was given and background ' \
                                             'range start must also be provided.'

        return issues

    def PyExec(self):
        from IndirectCommon import getInstrRun

        # Do setup
        self._setup()

        # Get input workspaces
        input_workspace_names = self._input_workspaces.getNames()

        # Lists of input and output workspaces
        q_workspaces = list()
        q2_workspaces = list()
        run_numbers = list()
        sample_param = list()

        progress = Progress(self, 0.0, 0.05, 3)
        q_ws, log_ws = None, None

        # Perform the ElasticWindow algorithms
        for input_ws in input_workspace_names:

            logger.information('Running ElasticWindow for workspace: %s' % input_ws)
            progress.report('ElasticWindow for workspace: %s' % input_ws)

            q_ws, q2_ws = ElasticWindow(InputWorkspace=input_ws,
                                        IntegrationRangeStart=self._integration_range_start,
                                        IntegrationRangeEnd=self._integration_range_end,
                                        OutputInQ="q_output", OutputInQSquared="q2_output",
                                        StoreInADS=False, EnableLogging=False)

            log_ws = Logarithm(InputWorkspace=q2_ws, OutputWorkspace="log_output",
                               StoreInADS=False, EnableLogging=False)

            q_workspaces.append(q_ws)
            q2_workspaces.append(log_ws)

            # Get the run number
            run_no = getInstrRun(input_ws)[1]
            run_numbers.append(run_no)

            # Get the sample environment unit
            sample, unit = self._get_sample_units(input_ws)
            if sample is not None:
                sample_param.append(sample)
            else:
                # No need to output a temperature workspace if there are no temperatures
                self._elt_ws_name = ''

        logger.information('Creating Q and Q^2 workspaces')
        progress.report('Creating Q workspaces')

        if len(input_workspace_names) == 1:
            # Just rename single workspaces
            q_workspace = q_ws
            q2_workspace = log_ws
        else:
            q_workspace = q_workspaces[0]
            q2_workspace = q2_workspaces[0]

            for idx in range(1, len(input_workspace_names)):
                q_workspace = _append_workspaces(q_workspace, q_workspaces[idx])
                q2_workspace = _append_workspaces(q2_workspace, q2_workspaces[idx])

        # Set the vertical axis units
        v_axis_is_sample = len(input_workspace_names) == len(sample_param)

        if v_axis_is_sample:
            logger.notice('Vertical axis is in units of %s' % unit)
            unit = (self._sample_log_name, unit)
        else:
            logger.notice('Vertical axis is in run number')
            unit = ('Run No', 'last 3 digits')

        # Create a new vertical axis for the Q and Q**2 workspaces
        q_ws_axis = NumericAxis.create(len(input_workspace_names))
        q_ws_axis.setUnit("Label").setLabel(unit[0], unit[1])

        q2_ws_axis = NumericAxis.create(len(input_workspace_names))
        q2_ws_axis.setUnit("Label").setLabel(unit[0], unit[1])

        # Set the vertical axis values
        for idx in range(0, len(input_workspace_names)):
            if v_axis_is_sample:
                q_ws_axis.setValue(idx, float(sample_param[idx]))
                q2_ws_axis.setValue(idx, float(sample_param[idx]))
            else:
                q_ws_axis.setValue(idx, float(run_numbers[idx][-3:]))
                q2_ws_axis.setValue(idx, float(run_numbers[idx][-3:]))

        # Add the new vertical axis to each workspace
        q_workspace.replaceAxis(1, q_ws_axis)
        q2_workspace.replaceAxis(1, q2_ws_axis)

        progress.report('Creating ELF workspaces')
        # Process the ELF workspace
        if self._elf_ws_name != '':
            logger.information('Creating ELF workspace')
            tranposed_q = _transpose_workspace(q_workspace)
            self.setProperty('OutputELF', _sort_workspace(tranposed_q))

        # Do temperature normalisation
        if self._elt_ws_name != '':
            logger.information('Creating ELT workspace')

            # If the ELF workspace was not created, create the ELT workspace
            # from the Q workspace. Else, clone the ELF workspace.
            if self._elf_ws_name == '':
                tranposed_q = _transpose_workspace(q_workspace)
                elt_workspace = _sort_workspace(tranposed_q)
            else:
                elt_workspace = _clone_workspace(self.getProperty('OutputELF').value)

            _normalize_by_index(elt_workspace, np.argmin(sample_param))
            self.setProperty('OutputELT', elt_workspace)

        # Set the output workspace
        self.setProperty('OutputInQ', q_workspace)
        self.setProperty('OutputInQSquared', q2_workspace)

    def _setup(self):
        """
        Gets algorithm properties.
        """

        self._sample_log_name = self.getPropertyValue('SampleEnvironmentLogName')
        self._sample_log_value = self.getPropertyValue('SampleEnvironmentLogValue')

        self._input_workspaces = self.getProperty('InputWorkspaces').value

        self._background_range_start = self.getProperty('BackgroundRangeStart').value
        self._background_range_end = self.getProperty('BackgroundRangeEnd').value

        if self._background_range_start != Property.EMPTY_DBL and self._background_range_end != Property.EMPTY_DBL:
            self._integration_range_start = self._background_range_start
            self._integration_range_end = self._background_range_end
        else:
            self._integration_range_start = self.getProperty('IntegrationRangeStart').value
            self._integration_range_end = self.getProperty('IntegrationRangeEnd').value


    def _get_sample_units(self, workspace):
        """
        Gets the sample environment units for a given workspace.

        @param workspace The workspace
        @returns sample in given units or None if not found
        """
        from IndirectCommon import getInstrRun

        instr, run_number = getInstrRun(workspace)

        facility = config.getFacility()
        pad_num = facility.instrument(instr).zeroPadding(int(run_number))
        zero_padding = '0' * (pad_num - len(run_number))

        run_name = instr + zero_padding + run_number
        log_filename = run_name.upper() + '.log'

        run = workspace.getRun()

        if self._sample_log_name == 'Position':
            # Look for sample changer position in logs in workspace
            if self._sample_log_name in run:
                tmp = run[self._sample_log_name].value
                value_action = {'last_value': lambda x: x[-1],
                                'average': lambda x: x.mean()}
                position = value_action['last_value'](tmp)
                if position == 0:
                    self._sample_log_name = 'Bot_Can_Top'
                if position == 1:
                    self._sample_log_name = 'Middle_Can_Top'
                if position == 2:
                    self._sample_log_name = 'Top_Can_Top'
            else:
                logger.information('Position not found in workspace.')

        if self._sample_log_name in run:
            # Look for sample unit in logs in workspace
            tmp = run[self._sample_log_name].value
            value_action = {'last_value': lambda x: x[-1],
                            'average': lambda x: x.mean()}
            sample = value_action[self._sample_log_value](tmp)
            unit = run[self._sample_log_name].units
            logger.information('%d %s found for run: %s' % (sample, unit, run_name))
            return sample, unit

        else:
            # Logs not in workspace, try loading from file
            logger.information('Log parameter not found in workspace. Searching for log file.')
            log_path = FileFinder.getFullPath(log_filename)

            if log_path != '':
                # Get temperature from log file
                LoadLog(Workspace=workspace, Filename=log_path)
                run_logs = workspace.getRun()
                if self._sample_log_name in run_logs:
                    tmp = run_logs[self._sample_log_name].value
                    sample = tmp[len(tmp) - 1]
                    unit = run[self._sample_log_name].units
                    logger.debug('%d %s found for run: %s' % (sample, unit, run_name))
                    return sample, unit
                else:
                    logger.warning('Log entry %s for run %s not found' % (self._sample_log_name, run_name))
            else:
                logger.warning('Log file for run %s not found' % run_name)

        # Can't find log file
        logger.warning('No sample units found for run: %s' % run_name)
        return None, None


# Register algorithm with Mantid
AlgorithmFactory.subscribe(ElasticWindowMultiple)
