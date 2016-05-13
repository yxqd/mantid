import numpy as np
import os
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from mantid.simpleapi import *


class PoldiPositionCalibration(object):
    def __init__(self, tableWorkspaceName):
        tableWorkspace = AnalysisDataService.retrieve(tableWorkspaceName)

        self.curves, self.x_min, self.x_max = self._get_curve_data(tableWorkspace)

    def _get_curve_data(self, calibrationTable):
        y_values = set(calibrationTable.column(2))
        x_values = set(calibrationTable.column(1))
        curves = dict([(x, dict([(y, {'raw': []}) for y in y_values])) for x in calibrationTable.column(0)])
        for i in range(calibrationTable.rowCount()):
            row = calibrationTable.row(i)

            curves[row['TwoTheta']][row['y0']]['raw'].append((row['x0'], row['a'], row['slope']))

            if len(curves[row['TwoTheta']][row['y0']]['raw']) == len(x_values):
                curves[row['TwoTheta']][row['y0']]['data'] = np.array(curves[row['TwoTheta']][row['y0']]['raw'])

        for two_theta, calibration_plot in curves.iteritems():
            for y0, data_dict in calibration_plot.iteritems():
                data = data_dict['data']

                d = {}

                poly_a = np.polyfit(data[:, 0], data[:, 1], 1, full=True)
                d['x0_a'] = (5.4311946 - poly_a[0][1]) / poly_a[0][0]

                poly_s = np.polyfit(data[:, 0], data[:, 2], 1, full=True)
                d['x0_slope'] = - poly_s[0][1] / poly_s[0][0]

                d['x0_diff'] = d['x0_a'] - d['x0_slope']
                d['x0_avg'] = (d['x0_a'] + d['x0_slope']) / 2

                data_dict['fit'] = d

        return curves, min(x_values), max(x_values)

    def generatePlots(self, path):
        if not os.path.exists(path):
            raise RuntimeError('The path you specified for saving does not exist: {0}'.format(path))

        for i, two_theta in enumerate(sorted(self.curves.keys())):
            figure = self._get_two_theta_plot(two_theta)
            figure.savefig(os.path.join(path, '{0}_{1}.png'.format(str(i), str(two_theta).replace('.', '_'))),
                           bbox_inches='tight')

    def generateCSV(self, filename):
        fh = open(filename, 'w')
        fh.write('#2theta,y0,x0_slope,x0_a,dx0,x0_avg\n')
        for two_theta in sorted(self.curves.keys()):
            y0_lines = self.curves[two_theta]
            for y0 in sorted(y0_lines.keys()):
                data = y0_lines[y0]
                fh.write(','.join((str(x) for x in (two_theta, y0,
                                                    data['fit']['x0_slope'], data['fit']['x0_a'],
                                                    data['fit']['x0_diff'], data['fit']['x0_avg']))) + '\n')

        fh.close()

    def _get_two_theta_plot(self, two_theta):
        axes_a, axes_slope, fig = self._get_formatted_two_theta_figure(two_theta)

        # Find minimum |x_0,slope - x_0,a| at which y0
        min_abs_diff = None
        x0_average = None
        y0_min = None

        for y0, data in self.curves[two_theta].iteritems():
            if min_abs_diff is None or abs(data['fit']['x0_diff']) < min_abs_diff:
                min_abs_diff = abs(data['fit']['x0_diff'])
                x0_average = data['fit']['x0_avg']
                y0_min = y0

            # Plot the y0-lines
            axes_a.plot(data['data'][:, 0], data['data'][:, 1], label=str(y0))
            axes_slope.plot(data['data'][:, 0], data['data'][:, 2], label=str(y0))

        x_margin = (self.x_max - self.x_min) * 0.05
        x_min_final = min(self.x_min, x0_average - x_margin)
        x_max_final = max(self.x_max, x0_average + x_margin)

        axes_a.axvline(x=x0_average, color='black')
        axes_a.set_xlim(xmin=x_min_final, xmax=x_max_final)
        axes_a.plot([x_min_final, x_max_final], [5.431196] * 2, ls='--', color='black')

        axes_slope.axvline(x=x0_average, ymin=-100000, ymax=100000, color='black')
        axes_slope.set_xlim(xmin=x_min_final, xmax=x_max_final)
        axes_slope.plot([x_min_final, x_max_final], [0] * 2, ls='--', color='black')
        axes_slope.legend()

        fig.text(0.5, 0.0,
                 '$y_0={0}, \\overline{{x}} = {1:.5f}, |\\Delta x_0| = {2:.5f}$'.format(y0_min, x0_average,
                                                                                        min_abs_diff), ha='center')

        fig.tight_layout()

        return fig

    def _get_formatted_two_theta_figure(self, two_theta):
        fig = Figure(figsize=(11.7, 8.3), dpi=90)
        canvas = FigureCanvas(fig)
        fig.suptitle('$2\\theta = {0}$'.format(str(two_theta)))
        axes_a = fig.add_subplot(121)
        axes_a.set_xlabel('$x_0$ in m')
        axes_a.set_ylabel('$a$ in $\AA{}$')
        axes_slope = fig.add_subplot(122)
        axes_slope.set_xlabel('$x_0$ in m')
        axes_slope.set_ylabel('Slope')
        return axes_a, axes_slope, fig
