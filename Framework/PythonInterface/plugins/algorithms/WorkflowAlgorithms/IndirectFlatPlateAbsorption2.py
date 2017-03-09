# pylint: disable=no-init,too-many-instance-attributes,too-many-branches
from mantid.simpleapi import *
from mantid.api import DataProcessorAlgorithm, AlgorithmFactory, PropertyMode, MatrixWorkspaceProperty, \
    WorkspaceGroupProperty, InstrumentValidator, WorkspaceUnitValidator, Progress
from mantid.kernel import StringListValidator, StringMandatoryValidator, IntBoundedValidator, \
    FloatBoundedValidator, Direction, logger, CompositeValidator


class IndirectFlatPlateAbsorption(DataProcessorAlgorithm):
    # Sample variables
    _sample_ws_name = None
    _sample_chemical_formula = None
    _sample_density_type = None
    _sample_density = 0.
    _sample_height = None
    _sample_width = None
    _sample_thickness = None

    # Container variables
    _can_ws_name = None
    _use_can_corrections = None
    _can_chemical_formula = None
    _can_density_type = None
    _can_density = 0.
    _can_front_thickness = None
    _can_back_thickness = None
    _can_scale = None

    _beam_height = None
    _beam_width = None
    _unit = None
    _emode = None
    _efixed = None
    _number_wavelengths = 10
    _events = 2000
    _abs_ws = None
    _ass_ws = None
    _acc_ws = None
    _output_ws = None

    def category(self):
        return "Workflow\\Inelastic;CorrectionFunctions\\AbsorptionCorrections;Workflow\\MIDAS"

    def summary(self):
        return "Calculates indirect absorption corrections for a flat sample shape."

    def version(self):
        return 2

    def PyInit(self):
        # Sample
        self.declareProperty(MatrixWorkspaceProperty('SampleWorkspace', '',
                                                     direction=Direction.Input),
                             doc='Sample workspace')
        self.declareProperty(name='SampleChemicalFormula', defaultValue='',
                             validator=StringMandatoryValidator(),
                             doc='Chemical formula for the sample')
        self.declareProperty(name='SampleDensityType', defaultValue='Mass Density',
                             validator=StringListValidator(['Mass Density', 'Number Density']),
                             doc='Sample density type.')
        self.declareProperty(name='SampleDensity', defaultValue=1.0,
                             validator=FloatBoundedValidator(0.0),
                             doc='Sample number density')
        self.declareProperty(name='SampleWidth', defaultValue=1.0,
                             validator=FloatBoundedValidator(0.0),
                             doc='Sample width')
        self.declareProperty(name='SampleThickness', defaultValue=0.5,
                             validator=FloatBoundedValidator(0.0),
                             doc='Sample thickness')
        self.declareProperty(name='SampleHeight', defaultValue=1.0,
                             validator=FloatBoundedValidator(0.0),
                             doc='Sample height')
        self.declareProperty(name='SampleAngle', defaultValue=0.0,
                             doc='Sample angle to beam')

        # Container
        self.declareProperty(MatrixWorkspaceProperty('ContainerWorkspace', '',
                                                     optional=PropertyMode.Optional,
                                                     direction=Direction.Input),
                             doc='Container workspace')
        self.declareProperty(name='UseContainerCorrections', defaultValue=False,
                             doc='Use Container corrections in subtraction')
        self.declareProperty(name='ContainerChemicalFormula', defaultValue='',
                             doc='Chemical formula for the Container')
        self.declareProperty(name='ContainerDensityType', defaultValue='Mass Density',
                             validator=StringListValidator(['Mass Density', 'Number Density']),
                             doc='Container density type.')
        self.declareProperty(name='ContainerDensity', defaultValue=1.0,
                             validator=FloatBoundedValidator(0.0),
                             doc='Container number density')
        self.declareProperty(name='ContainerFrontThickness', defaultValue=0.1,
                             validator=FloatBoundedValidator(0.0),
                             doc='Container front thickness')
        self.declareProperty(name='ContainerBackThickness', defaultValue=0.1,
                             validator=FloatBoundedValidator(0.0),
                             doc='Container back thickness')
        self.declareProperty(name='ContainerScaleFactor', defaultValue=1.0,
                             validator=FloatBoundedValidator(0.0),
                             doc='Scale factor to multiply Container data')

        # Beam size
        self.declareProperty(name='BeamHeight', defaultValue=1.0,
                             validator=FloatBoundedValidator(0.0),
                             doc='Height of the beam (cm)')
        self.declareProperty(name='BeamWidth', defaultValue=1.0,
                             validator=FloatBoundedValidator(0.0),
                             doc='Width of the beam (cm)')

        # Monte Carlo
        self.declareProperty(name='NumberOfWavelengthPoints', defaultValue=10,
                             validator=IntBoundedValidator(1),
                             doc='Number of wavelengths for calculation')
        self.declareProperty(name='EventsPerPoint', defaultValue=1000,
                             validator=IntBoundedValidator(0),
                             doc='Number of neutron events')

        # Output
        self.declareProperty(MatrixWorkspaceProperty('OutputWorkspace', '',
                                                     direction=Direction.Output),
                             doc='The output corrected workspace')

        self.declareProperty(WorkspaceGroupProperty('CorrectionsWorkspace', '',
                                                    direction=Direction.Output,
                                                    optional=PropertyMode.Optional),
                             doc='The workspace group to save correction factors')

    def PyExec(self):

        # Set up progress reporting
        n_prog_reports = 2
        if self._can_ws_name is not None:
            n_prog_reports += 1
        prog = Progress(self, 0.0, 1.0, n_prog_reports)
        delete_alg = self.createChildAlgorithm("DeleteWorkspace", enableLogging=False)
        multiply_alg = self.createChildAlgorithm("Multiply", enableLogging=False)
        divide_alg = self.createChildAlgorithm("Divide", enableLogging=False)
        minus_alg = self.createChildAlgorithm("Minus", enableLogging=False)
        convert_unit_alg = self.createChildAlgorithm("ConvertUnits", enableLogging=False)
        clone_alg = self.createChildAlgorithm("CloneWorkspace", enableLogging=False)
        group_alg = self.createChildAlgorithm("GroupWorkspaces", enableLogging=False)

        sample_wave_ws = '__sam_wave'
        if self._unit == 'Wavelength':
            clone_alg.setProperty("InputWorkspace", self._sample_ws_name)
            clone_alg.setProperty("OutputWorkspace", sample_wave_ws)
            clone_alg.execute()
            mtd.addOrReplace(sample_wave_ws, clone_alg.getProperty("OutputWorkspace").value)
        else:
            convert_unit_alg.setProperty("InputWorkspace", self._sample_ws_name)
            convert_unit_alg.setProperty("OutputWorkspace", sample_wave_ws)
            convert_unit_alg.setProperty("Target", 'Wavelength')
            convert_unit_alg.setProperty("EMode", self._emode)
            if self._emode == 'Indirect':
                convert_unit_alg.setProperty("EFixed", self._efixed)
            convert_unit_alg.execute()
            mtd.addOrReplace(sample_wave_ws, convert_unit_alg.getProperty("OutputWorkspace").value)

        prog.report('Calculating sample corrections')
        FlatPlateMonteCarloAbsorption(InputWorkspace=sample_wave_ws,
                                      OutputWorkspace=self._ass_ws,
                                      ChemicalFormula=self._sample_chemical_formula,
                                      DensityType=self._sample_density_type,
                                      Density=self._sample_density,
                                      Height=self._sample_height,
                                      Width=self._sample_width,
                                      Thickness=self._sample_thickness,
                                      Angle=self._sample_angle,
                                      Center=0.,
                                      BeamHeight=self._beam_height,
                                      BeamWidth=self._beam_width,
                                      EventsPerPoint=self._events,
                                      NumberOfWavelengthPoints=self._number_wavelengths,
                                      Interpolation=self._interpolation)
        group = self._ass_ws

        if self._can_ws_name is not None:
            can1_wave_ws = '__can1_wave'
            can2_wave_ws = '__can2_wave'
            convert_unit_alg.setProperty("InputWorkspace", self._can_ws_name)
            convert_unit_alg.setProperty("OutputWorkspace", can1_wave_ws)
            convert_unit_alg.setProperty("Target", 'Wavelength')
            convert_unit_alg.setProperty("EMode", self._emode)
            convert_unit_alg.setProperty("EFixed", self._efixed)
            convert_unit_alg.execute()
            mtd.addOrReplace(can1_wave_ws, convert_unit_alg.getProperty("OutputWorkspace").value)

            if self._can_scale != 1.0:
                logger.information('Scaling container by: %s' % self._can_scale)
                scale_alg = self.createChildAlgorithm("Scale", enableLogging=False)
                scale_alg.setProperty("InputWorkspace", can1_wave_ws)
                scale_alg.setProperty("OutputWorkspace", can1_wave_ws)
                scale_alg.setProperty("Factor", self._can_scale)
                scale_alg.setProperty("Operation", 'Multiply')
                scale_alg.execute()
            clone_alg = self.createChildAlgorithm("CloneWorkspace", enableLogging=False)
            clone_alg.setProperty("InputWorkspace", can1_wave_ws)
            clone_alg.setProperty("OutputWorkspace", can2_wave_ws)
            clone_alg.execute()
            mtd.addOrReplace(can2_wave_ws, clone_alg.getProperty("OutputWorkspace").value)

            if self._use_can_corrections:
                prog.report('Calculating container corrections')
                divide_alg.setProperty("LHSWorkspace", sample_wave_ws)
                divide_alg.setProperty("RHSWorkspace", self._ass_ws)
                divide_alg.setProperty("OutputWorkspace", sample_wave_ws)
                divide_alg.execute()

                offset_front = 0.5 * (float(self._can_front_thickness) + float(self._sample_thickness))
                FlatPlateMonteCarloAbsorption(InputWorkspace=can1_wave_ws,
                                              OutputWorkspace='__Acc1',
                                              ChemicalFormula=self._can_chemical_formula,
                                              DensityType=self._can_density_type,
                                              Density=self._can_density,
                                              Height=self._sample_height,
                                              Width=self._sample_width,
                                              Thickness=self._can_front_thickness,
                                              Center=-offset_front,
                                              Angle=self._sample_angle,
                                              BeamHeight=self._beam_height,
                                              BeamWidth=self._beam_width,
                                              EventsPerPoint=self._events,
                                              NumberOfWavelengthPoints=self._number_wavelengths,
                                              Interpolation=self._interpolation)

                offset_back = 0.5 * (float(self._can_back_thickness) + float(self._sample_thickness))
                FlatPlateMonteCarloAbsorption(InputWorkspace=can2_wave_ws,
                                              OutputWorkspace='__Acc2',
                                              ChemicalFormula=self._can_chemical_formula,
                                              DensityType=self._can_density_type,
                                              Density=self._can_density,
                                              Height=self._sample_height,
                                              Width=self._sample_width,
                                              Thickness=self._can_back_thickness,
                                              Center=offset_back,
                                              Angle=self._sample_angle,
                                              BeamHeight=self._beam_height,
                                              BeamWidth=self._beam_width,
                                              EventsPerPoint=self._events,
                                              NumberOfWavelengthPoints=self._number_wavelengths,
                                              Interpolation=self._interpolation)

                multiply_alg.setProperty("LHSWorkspace", '__Acc1')
                multiply_alg.setProperty("RHSWorkspace", '__Acc2')
                multiply_alg.setProperty("OutputWorkspace", self._acc_ws)
                multiply_alg.execute()
                mtd.addOrReplace(self._acc_ws, multiply_alg.getProperty("OutputWorkspace").value)
                delete_alg.setProperty("Workspace", '__Acc1')
                delete_alg.execute()
                delete_alg.setProperty("Workspace", '__Acc2')
                delete_alg.execute()

                divide_alg.setProperty("LHSWorkspace", can1_wave_ws)
                divide_alg.setProperty("RHSWorkspace", self._acc_ws)
                divide_alg.setProperty("OutputWorkspace", can1_wave_ws)
                divide_alg.execute()
                minus_alg.setProperty("LHSWorkspace", sample_wave_ws)
                minus_alg.setProperty("RHSWorkspace", can1_wave_ws)
                minus_alg.setProperty("OutputWorkspace", sample_wave_ws)
                minus_alg.execute()
                group += ',' + self._acc_ws

            else:
                prog.report('Calculating container scaling')
                minus_alg.setProperty("LHSWorkspace", sample_wave_ws)
                minus_alg.setProperty("RHSWorkspace", can1_wave_ws)
                minus_alg.setProperty("OutputWorkspace", sample_wave_ws)
                minus_alg.execute()
                divide_alg.setProperty("LHSWorkspace", sample_wave_ws)
                divide_alg.setProperty("RHSWorkspace", self._ass_ws)
                divide_alg.setProperty("OutputWorkspace", sample_wave_ws)
                divide_alg.execute()

            delete_alg.setProperty("Workspace", can1_wave_ws)
            delete_alg.execute()
            delete_alg.setProperty("Workspace", can2_wave_ws)
            delete_alg.execute()

        else:
            divide_alg.setProperty("LHSWorkspace", sample_wave_ws)
            divide_alg.setProperty("RHSWorkspace", self._ass_ws)
            divide_alg.setProperty("OutputWorkspace", sample_wave_ws)
            divide_alg.execute()

        if self._unit == 'Wavelength':
            clone_alg.setProperty("InputWorkspace", sample_wave_ws)
            clone_alg.setProperty("OutputWorkspace", self._output_ws)
            clone_alg.execute()
            mtd.addOrReplace(self._output_ws, clone_alg.getProperty("OutputWorkspace").value)
        else:
            convert_unit_alg.setProperty("InputWorkspace", sample_wave_ws)
            convert_unit_alg.setProperty("OutputWorkspace", self._output_ws)
            convert_unit_alg.setProperty("Target", self._unit)
            convert_unit_alg.setProperty("EMode", self._emode)
            if self._emode == 'Inelastic':
                convert_unit_alg.setProperty("EFixed", self._efixed)
            convert_unit_alg.execute()
            mtd.addOrReplace(self._output_ws, convert_unit_alg.getProperty("OutputWorkspace").value)
        delete_alg.setProperty("Workspace", sample_wave_ws)
        delete_alg.execute()

        prog.report('Recording sample logs')
        sample_log_workspaces = [self._output_ws, self._ass_ws]
        sample_logs = [('sample_shape', 'flatplate'),
                       ('sample_filename', self._sample_ws_name),
                       ('sample_height', self._sample_height),
                       ('sample_width', self._sample_width),
                       ('sample_thickness', self._sample_thickness)]

        if self._can_ws_name is not None:
            sample_logs.append(('container_filename', self._can_ws_name))
            sample_logs.append(('container_scale', self._can_scale))
            if self._use_can_corrections:
                sample_log_workspaces.append(self._acc_ws)
                sample_logs.append(('container_front_thickness', self._can_front_thickness))
                sample_logs.append(('container_back_thickness', self._can_back_thickness))

        log_names = [item[0] for item in sample_logs]
        log_values = [item[1] for item in sample_logs]

        add_sample_log_alg = self.createChildAlgorithm("AddSampleLogMultiple", enableLogging=False)
        for ws_name in sample_log_workspaces:
            add_sample_log_alg.setProperty("Workspace", ws_name)
            add_sample_log_alg.setProperty("LogNames", log_names)
            add_sample_log_alg.setProperty("LogValues", log_values)
            add_sample_log_alg.execute()

        self.setProperty('OutputWorkspace', self._output_ws)

        # Output the Ass workspace if it is wanted, delete if not
        if self._abs_ws == '':
            delete_alg.setProperty("Workspace", self._ass_ws)
            delete_alg.execute()
            if self._can_ws_name is not None and self._use_can_corrections:
                delete_alg.setProperty("Workspace", self._acc_ws)
                delete_alg.execute()
        else:
            group_alg.setProperty("InputWorkspaces", group)
            group_alg.setProperty("OutputWorkspace", self._abs_ws)
            group_alg.execute()
            mtd.addOrReplace(self._abs_ws, group_alg.getProperty("OutputWorkspace").value)
            self.setProperty('CorrectionsWorkspace', self._abs_ws)

    def _setup(self):
        """
        Get algorithm properties.
        """

        self._sample_ws_name = self.getPropertyValue('SampleWorkspace')
        ws = mtd[self._sample_ws_name]
        axis = ws.getAxis(0)
        self._unit = axis.getUnit().unitID()
        logger.information('Input X-unit is %s' % self._unit)
        if self._unit == 'dSpacing':
            self._emode = 'Elastic'
        else:
            self._emode = str(ws.getEMode())
        logger.information('Input Emode is %s' % self._emode)
        if self._emode == 'Indirect':
            self._efixed = self._get_Efixed()
            logger.information('Input Efixed is %f' % self._efixed)

        self._beam_height = self.getProperty('BeamHeight').value
        self._beam_width = self.getProperty('BeamWidth').value

        self._sample_chemical_formula = self.getPropertyValue('SampleChemicalFormula')
        self._sample_density_type = self.getProperty('SampleDensityType').value
        self._sample_density = self.getProperty('SampleDensity').value
        sample_height = self.getProperty('SampleHeight').value
        if sample_height == '':
            self._sample_height = self._beam_height
        else:
            self._sample_height = sample_height
        sample_width = self.getProperty('SampleWidth').value
        if sample_width == '':
            self._sample_width = self._beam_width
        else:
            self._sample_width = sample_width
        self._sample_thickness = self.getProperty('SampleThickness').value
        self._sample_angle = self.getProperty('SampleAngle').value

        self._can_ws_name = self.getPropertyValue('ContainerWorkspace')
        if self._can_ws_name == '':
            self._can_ws_name = None
        self._use_can_corrections = self.getProperty('UseContainerCorrections').value
        self._can_chemical_formula = self.getPropertyValue('ContainerChemicalFormula')
        self._can_density_type = self.getProperty('ContainerDensityType').value
        self._can_density = self.getProperty('ContainerDensity').value
        self._can_front_thickness = self.getProperty('ContainerFrontThickness').value
        self._can_back_thickness = self.getProperty('ContainerBackThickness').value
        self._can_scale = self.getProperty('ContainerScaleFactor').value

        self._number_wavelengths = self.getProperty('NumberOfWavelengthPoints').value
        self._events = self.getProperty('EventsPerPoint').value
        self._interpolation = 'CSpline'

        self._output_ws = self.getPropertyValue('OutputWorkspace')

        self._abs_ws = self.getPropertyValue('CorrectionsWorkspace')
        if self._abs_ws == '':
            self._ass_ws = '__ass'
            self._acc_ws = '__acc'
        else:
            self._ass_ws = self._abs_ws + '_ass'
            self._acc_ws = self._abs_ws + '_acc'

    def validateInputs(self):
        """
        Validate algorithm options.
        """

        self._setup()
        issues = dict()

        if self._use_can_corrections and self._can_chemical_formula == '':
            issues['CanChemicalFormula'] = 'Must be set to use can corrections'

        if self._use_can_corrections and self._can_ws_name is None:
            issues['UseCanCorrections'] = 'Must specify a can workspace to use can corrections'

        return issues

    def _get_Emode(self):
        inst = mtd[self._sample_ws_name].getInstrument()

        if inst.hasParameter('Emode'):
            return inst.getStringParameter('Emode')[0]
        raise ValueError('No Emode parameter found')

    def _get_Efixed(self):
        inst = mtd[self._sample_ws_name].getInstrument()

        if inst.hasParameter('Efixed'):
            return inst.getNumberParameter('EFixed')[0]

        if inst.hasParameter('analyser'):
            analyser_name = inst.getStringParameter('analyser')[0]
            analyser_comp = inst.getComponentByName(analyser_name)

            if analyser_comp is not None and analyser_comp.hasParameter('Efixed'):
                return analyser_comp.getNumberParameter('EFixed')[0]

        raise ValueError('No Efixed parameter found')


# Register algorithm with Mantid
AlgorithmFactory.subscribe(IndirectFlatPlateAbsorption)
