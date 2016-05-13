
.. algorithm::

.. summary::

.. alias::

.. properties::

Description
-----------

This algorithm is used for the calibration of the POLDI instrument parameters. There are two parameters related to the
chopper, :math:`t_0` and :math:`t_{const}`, and three parameters related to the detector position, the center of the
detector circle (:math:`x_0` and :math:`y_0`) and the angle between the primary beam and the central detector element
(:math:`2\theta_0`).

The calibration is performed using a silicon standard with :math:`a = 5.4311946 \AA{}`. For a full calibration of all
five parameters, the standard must be measured at three different chopper speeds. First, the chopper speed related
parameters are determined. An offset in :math:`t_{const}` leads to a :math:`Q`-dependence of the observed lattice
parameter, while an offset in :math:`t_0` leads to a different :math:`Q`-dependence for each chopper speed.
These two parameters are optimized independently from each other, after this the lattice parameters determined from each
single peak of the silicon standard should be the same, altough they may deviate from the absolute value.

Calibration of the detector position is performed using data collected at the highest chopper speed. For this, a special
function is used, which has a "slope" parameter that is very sensitive to tilts in the detector. The parameter must as
close to 0 as possible. In addition, the absolute value of the lattice parameter must be as close to the certified
standard parameter as possible. Currently the optimum parameters are found by fitting the slope- and lattice parameter
for all parameter combination is a certain range. At each :math:`2\theta`-value, slope and lattice parameter are plotted
against :math:`x_0`, with one line for each of the :math:`y_0`-values. To find the correct parameter combination one has
to find the combination of :math:`2\theta` and :math:`y_0` where the intersections of the plotted lines with
:math:`y=0` for the slope and :math:`y=5.43211946` for :math:`a` are as close as possible.

Usually more than one iteration of this process has to be performed, so that new timing parameters are found with the
new position parameters and so on.

Usage
-----

.. note::

    This algorithm is intended to be used by the POLDI beamline scientists. The examples here are only guides how to
    use the algorithm, actually performing the full calibration is very time consuming and only possible with the
    calibration data available.

First of all, the calibration data for three different speed have to be loaded using :ref:`algm-PoldiLoadRuns`, so that
they have a valid instrument (the latest valid instrument definition and parameters). Furthermore, the peaks for Si are
required for indexing:


.. code-block::

    PoldiLoadRuns(2015, 6495, 6499, 5, OutputWorkspace='calibration')
    PoldiLoadRuns(2015, 6500, 6505, 6, OutputWorkspace='calibration')
    PoldiLoadRuns(2015, 6506, 6517, 12, OutputWorkspace='calibration')
    PoldiCreatePeaksFromCell(SpaceGroup='F d -3 m',
                             Atoms='Si 0 0 0 1.0 0.01',
                             a=5.4311946, LatticeSpacingMin=0.7, OutputWorkspace='Si')

Next, the timing parameters have to be determined using the respective mode:

.. code-block::

    PoldiCalibration(InputWorkspace='calibration', ExpectedPeaks='Si', CalibrationMode='Timing',
                     OutputWorkspace='calibration_timing_1')

The result of this is a TableWorkspace containing the two calibrated timing related parameters. These have to be used
in the next iteration which is a calibration step for the detector position. In addition, the ranges for the
three position parameters have to be specified in the format "start, stop, number of steps" in the order
:math:`2\theta`, :math:`x_0` (in m), :math:`y_0`. The number of combinations grows quickly, so it makes sense to
start with a few values first and make the ranges smaller afterwards when it's more clear where the parameters are
approximately.

.. code-block::

    PoldiCalibration(InputWorkspace='calibration', ExpectedPeaks='Si', CalibrationMode='Position',
                     InitialParameters='-0.037593, -11.141',
                     ParameterRanges='90.0,91.5,3; -0.900,-0.800,3; -0.900,-0.800,3'
                     OutputWorkspace='calibration_position_1')

The output of this algorithm is a table workspace that contains the slope and lattice parameter (along with their
errors) for each parameter combination. Furthermore it contains the average relative FWHM of the peaks, which should be
as small as possible. Usually, these values should be in the order of 0.002.




.. categories::

.. sourcelink::

