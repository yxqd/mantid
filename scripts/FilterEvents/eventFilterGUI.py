#pylint: disable=invalid-name, too-many-lines, too-many-instance-attributes
import numpy

from ui_MainWindow import Ui_MainWindow #import line for the UI python class
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from matplotlib.pyplot import setp

import mantid
import mantid.simpleapi as api
import mantid.kernel
from mantid.kernel import Logger
from mantid.simpleapi import AnalysisDataService

from mantid.kernel import ConfigService

import os

HUGE_FAST = 10000
HUGE_PARALLEL = 100000
MAXTIMEBINSIZE = 3000

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s


class MyPopErrorMsg(QWidget):
    """ Pop up dialog window
    """

    def __init__(self):
        """ Init
        """
        import ui_ErrorMessage as errui
        QWidget.__init__(self)

        self.ui = errui.Ui_Dialog()
        self.ui.setupUi(self)

        QtCore.QObject.connect(self.ui.pushButton_quit, QtCore.SIGNAL('clicked()'), self.quit)

    def setMessage(self, errmsg):
        """ Set message
        """
        self.ui.label_errmsg.setWordWrap(True)
        self.ui.label_errmsg.setText(errmsg)

        return

    def quit(self):
        """ Quit
        """
        self.close()

        return

    def XpaintEvent(self, _):
        """ ???
        """
        import ui_ErrorMessage as errui

        self.ui = errui.Ui_Dialog()
        self.ui.setupUi(self)

        return


class MainWindow(QtGui.QMainWindow):
    """ Class of Main Window (top)
    """

    _errMsgWindow = None

    def __init__(self, parent=None):
        """ Intialization and set up
        """
        # Base class
        QtGui.QMainWindow.__init__(self,parent)

        # Mantid configuration
        config = ConfigService.Instance()
        self._instrument = config["default.instrument"]

        # Central widget
        self.centralwidget = QtGui.QWidget(self)

        # UI Window (from Qt Designer)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.mainplot = self.ui.graphicsView.getPlot()

        # Do initialize plotting
        vecx, vecy, xlim, ylim = self.computeMock()

        self.mainline = self.ui.mainplot.plot(vecx, vecy, 'r-')

        leftx = [xlim[0], xlim[0]]
        lefty = [ylim[0], ylim[1]]
        self.leftslideline = self.ui.mainplot.plot(leftx, lefty, 'b--')
        rightx = [xlim[1], xlim[1]]
        righty = [ylim[0], ylim[1]]
        self.rightslideline = self.ui.mainplot.plot(rightx, righty, 'g--')
        upperx = [xlim[0], xlim[1]]
        uppery = [ylim[1], ylim[1]]
        self.upperslideline = self.ui.mainplot.plot(upperx, uppery, 'b--')
        lowerx = [xlim[0], xlim[1]]
        lowery = [ylim[0], ylim[0]]
        self.lowerslideline = self.ui.mainplot.plot(lowerx, lowery, 'g--')

        self.ui.graphicsView.mpl_connect('button_press_event', self.on_mouseDownEvent)

        # Set up horizontal slide (integer) and string value
        self._leftSlideValue = 0
        self._rightSlideValue = 99

        self.ui.horizontalSlider.setRange(0, 100)
        self.ui.horizontalSlider.setValue(self._leftSlideValue)
        self.ui.horizontalSlider.setTracking(True)
        self.ui.horizontalSlider.setTickPosition(QSlider.NoTicks)
        self.connect(self.ui.horizontalSlider, SIGNAL('valueChanged(int)'), self.move_leftSlider)

        self.ui.horizontalSlider_2.setRange(0, 100)
        self.ui.horizontalSlider_2.setValue(self._rightSlideValue)
        self.ui.horizontalSlider_2.setTracking(True)
        self.ui.horizontalSlider_2.setTickPosition(QSlider.NoTicks)
        self.connect(self.ui.horizontalSlider_2, SIGNAL('valueChanged(int)'), self.move_rightSlider)

        # self.connect(self.ui.lineEdit_3, QtCore.SIGNAL("textChanged(QString)"),
        #         self.set_startTime)
        self.ui.lineEdit_startTime.setValidator(QtGui.QDoubleValidator(self.ui.lineEdit_startTime))
        self.connect(self.ui.pushButton_setT0, QtCore.SIGNAL("clicked()"), self.set_startTime)
        # self.connect(self.ui.lineEdit_4, QtCore.SIGNAL("textChanged(QString)"),
        #         self.set_stopTime)
        self.ui.lineEdit_stopTime.setValidator(QtGui.QDoubleValidator(self.ui.lineEdit_stopTime))
        self.connect(self.ui.pushButton_setTf, QtCore.SIGNAL("clicked()"), self.set_stopTime)

        # File loader
        self.scan_event_worksapces()
        self.connect(self.ui.pushButton_refreshWS, SIGNAL('clicked()'), self.scan_event_worksapces)
        self.connect(self.ui.pushButton_browse, SIGNAL('clicked()'), self.browse_File)
        self.connect(self.ui.pushButton_load, SIGNAL('clicked()'), self.load_file)
        self.connect(self.ui.pushButton_3, SIGNAL('clicked()'), self.use_existWS)

        # Set up time
        # self.ui.lineEdit_3.setValidator(QtGui.QDoubleValidator(self.ui.lineEdit_3))
        # self.ui.lineEdit_4.setValidator(QtGui.QDoubleValidator(self.ui.lineEdit_4))

        # Filter by time
        self.connect(self.ui.pushButton_filterTime, SIGNAL('clicked()'), self.filter_by_time)

        # Filter by log value
        self.ui.lineEdit_minValue.setValidator(QtGui.QDoubleValidator(self.ui.lineEdit_minValue))
        self.ui.lineEdit_maxValue.setValidator(QtGui.QDoubleValidator(self.ui.lineEdit_maxValue))
        self.connect(self.ui.lineEdit_minValue, QtCore.SIGNAL("textChanged(QString)"),
                     self.set_min_log_value)
        self.connect(self.ui.lineEdit_maxValue, QtCore.SIGNAL("textChanged(QString)"),
                     self.set_max_log_value)


        self.ui.lineEdit_7.setValidator(QtGui.QDoubleValidator(self.ui.lineEdit_7))
        self.ui.lineEdit_8.setValidator(QtGui.QDoubleValidator(self.ui.lineEdit_8))
        self.ui.lineEdit_9.setValidator(QtGui.QDoubleValidator(self.ui.lineEdit_9))



        dirchangeops = ["Both", "Increase", "Decrease"]
        self.ui.comboBox_4.addItems(dirchangeops)

        logboundops = ["Centre", "Left"]
        self.ui.comboBox_5.addItems(logboundops)

        self.connect(self.ui.pushButton_4, SIGNAL('clicked()'), self.plotLogValue)

        self.connect(self.ui.pushButton_filterLog, SIGNAL('clicked()'), self.filter_by_log_value)

        #Set up help button
        self.connect(self.ui.helpBtn, QtCore.SIGNAL('clicked()'), self.helpClicked)

        # Set up vertical slide
        self._upperSlideValue = 99
        self._lowerSlideValue = 0

        self.ui.verticalSlider.setRange(0, 100)
        self.ui.verticalSlider.setValue(self._upperSlideValue)
        self.ui.verticalSlider.setTracking(True)
        self.connect(self.ui.verticalSlider, SIGNAL('valueChanged(int)'), self.move_upper_slider)

        self.ui.verticalSlider_2.setRange(0, 100)
        self.ui.verticalSlider_2.setValue(self._lowerSlideValue)
        self.ui.verticalSlider_2.setTracking(True)
        self.connect(self.ui.verticalSlider_2, SIGNAL('valueChanged(int)'), self.move_lowerSlider)

        # Set up for filtering (advanced setup)
        self._tofcorrection = False
        self.ui.checkBox_fastLog.setChecked(False)
        self.ui.checkBox_filterByPulse.setChecked(False)
        self.ui.checkBox_from1.setChecked(False)
        self.ui.checkBox_groupWS.setChecked(True)

        self.connect(self.ui.comboBox_tofCorr, SIGNAL('currentIndexChanged(int)'), self.showHideEi)
        self.connect(self.ui.pushButton_refreshCorrWSList, SIGNAL('clicked()'),  self._searchTableWorkspaces)

        self.ui.lineEdit_Ei.setValidator(QtGui.QDoubleValidator(self.ui.lineEdit_Ei))

        self.ui.label_Ei.hide()
        self.ui.lineEdit_Ei.hide()
        self.ui.label_Ei_2.hide()
        self.ui.comboBox_corrWS.hide()
        self.ui.pushButton_refreshCorrWSList.hide()

        # Error message
        # self.connect(self.ui.pushButton_clearerror, SIGNAL('clicked()'), self._clearErrorMsg)
        # self.ui.plainTextEdit_ErrorMsg.setReadOnly(True)
        # self.ui.label_error.hide()

        # Set up for workspaces
        self._dataWS = None
        self._sampleLogNames = []
        self._sampleLog = None

        # Side information
        self.ui.label_mean.hide()
        self.ui.label_meanvalue.hide()
        self.ui.label_avg.hide()
        self.ui.label_timeAvgValue.hide()
        self.ui.label_freq.hide()
        self.ui.label_freqValue.hide()
        self.ui.label_logname.hide()
        self.ui.label_lognamevalue.hide()
        self.ui.label_logsize.hide()
        self.ui.label_logsizevalue.hide()

        # Default
        self._defaultdir = os.getcwd()

        # self.ui.InputVal.setValidator(QtGui.QDoubleValidator(self.ui.InputVal))

        # QtCore.QObject.connect(self.ui.convert, QtCore.SIGNAL("clicked()"), self.convert )
        # QtCore.QObject.connect(self.ui.inputUnits, QtCore.SIGNAL("currentIndexChanged(QString)"), self.setInstrumentInputs )
        # QtCore.QObject.connect(self.ui.outputUnits, QtCore.SIGNAL("currentIndexChanged(QString)"), self.setInstrumentInputs )
        # self.setInstrumentInputs()

        ##defaults

        #register startup
        mantid.UsageService.registerFeatureUsage("Interface","EventFilter",False)

        return

    def on_mouseDownEvent(self, event):
        """ Respond to pick up a value with mouse down event
        """
        x = event.xdata
        y = event.ydata

        if x is not None and y is not None:
            msg = "You've clicked on a bar with coords:\n %f, %f" % (x, y)
            QMessageBox.information(self, "Click!", msg)

        return

    def computeMock(self):
        """ Compute vecx and vecy as mocking
        """
        x0 = 0.
        xf = 1.
        dx = 0.1

        vecx = []
        vecy = []

        x = x0
        while x < xf:
            y = 0.0
            vecx.append(x)
            vecy.append(y)
            x += dx

        xlim = [x0, xf]
        ylim = [-1., 1]

        return (vecx, vecy, xlim, ylim)

    def move_leftSlider(self):
        """ Re-setup left range line in figure.
        Triggered by a change in Qt Widget.  NO EVENT is required.
        """
        newx = self.ui.horizontalSlider.value()
        if newx <= self._rightSlideValue and newx != self._leftSlideValue:
            # Allowed value: move the value bar
            self._leftSlideValue = newx

            # Move the vertical line
            xlim = self.ui.mainplot.get_xlim()
            newx = xlim[0] + newx*(xlim[1] - xlim[0])*0.01
            leftx = [newx, newx]
            lefty = self.ui.mainplot.get_ylim()
            setp(self.leftslideline, xdata=leftx, ydata=lefty)

            self.ui.graphicsView.draw()

            # Change value
            self.ui.lineEdit_3.setText(str(newx))

        else:
            # Reset the value to original value
            self.ui.horizontalSlider.setValue(self._leftSlideValue)

        return

    def set_startTime(self):
        """ Set the starting time and left slide bar
        """
        inps = str(self.ui.lineEdit_3.text())
        info_msg = "Starting time = %s" % (inps)
        Logger("Filter_Events").information(info_msg)

        xlim = self.ui.mainplot.get_xlim()
        if inps == "":
            # Empty. Use default
            newtime0 = xlim[0]
        else:
            newtime0 = float(inps)

        # Convert to integer slide value
        ileftvalue = int( (newtime0-xlim[0])/(xlim[1] - xlim[0])*100 )
        debug_msg = "iLeftSlide = %s" % str(ileftvalue)
        Logger("Filter_Events").debug(debug_msg)

        # Skip if same as origina
        if ileftvalue == self._leftSlideValue:
            return

        # Set the value if out of range
        resetT = True
        if ileftvalue < 0:
            # Minimum value as 0
            ileftvalue = 0
        elif ileftvalue > self._rightSlideValue:
            # Maximum value as right slide value
            ileftvalue = self._rightSlideValue
        else:
            resetT = False

        if resetT is True:
            newtime0 = xlim[0] + ileftvalue*(xlim[1]-xlim[0])*0.01
        info_msg = "Corrected iLeftSlide = %s (vs. right = %s)" % (str(ileftvalue), str(self._rightSlideValue))
        Logger("Filter_Events").information(info_msg)

        # Move the slide bar (left)
        self._leftSlideValue = ileftvalue

        # Move the vertical line
        leftx = [newtime0, newtime0]
        lefty = self.ui.mainplot.get_ylim()
        setp(self.leftslideline, xdata=leftx, ydata=lefty)

        self.ui.graphicsView.draw()

        # Set the value to left slider
        self.ui.horizontalSlider.setValue(self._leftSlideValue)
        # Reset the value of line edit
        if resetT is True:
            self.ui.lineEdit_3.setText(str(newtime0))

        return

    def move_rightSlider(self):
        """ Re-setup left range line in figure.
        Triggered by a change in Qt Widget.  NO EVENT is required.
        """
        newx = self.ui.horizontalSlider_2.value()
        if newx >= self._leftSlideValue and newx != self._rightSlideValue:
            # Allowed value: move the value bar
            self._rightSlideValue = newx

            xlim = self.ui.mainplot.get_xlim()
            newx = xlim[0] + newx*(xlim[1] - xlim[0])*0.01
            leftx = [newx, newx]
            lefty = self.ui.mainplot.get_ylim()
            setp(self.rightslideline, xdata=leftx, ydata=lefty)

            self.ui.graphicsView.draw()

            # Change value
            self.ui.lineEdit_4.setText(str(newx))

        else:
            # Reset the value
            self.ui.horizontalSlider_2.setValue(self._rightSlideValue)

        return

    def set_stopTime(self):
        """ Set the starting time and left slide bar
        """
        inps = str(self.ui.lineEdit_4.text())
        info_msg = "Stopping time = %s" % (inps)
        Logger("Filter_Events").information(info_msg)

        xlim = self.ui.mainplot.get_xlim()
        if inps == "":
            # Empty. Use default
            newtimef = xlim[1]
        else:
            # Parse
            newtimef = float(inps)

        # Convert to integer slide value
        irightvalue = int( (newtimef-xlim[0])/(xlim[1] - xlim[0])*100 )
        info_msg = "iRightSlide = %s" % str(irightvalue)
        Logger("Filter_Events").information(info_msg)

        # Return if no change
        if irightvalue == self._rightSlideValue:
            return

        # Correct value
        resetT = True
        if irightvalue > 100:
            irightvalue = 100
        elif irightvalue < self._leftSlideValue:
            irightvalue = self._leftSlideValue
        else:
            resetT = False

        if resetT is True:
            newtimef = xlim[0] + irightvalue*(xlim[1]-xlim[0])*0.01

        # Move the slide bar (right)
        self._rightSlideValue = irightvalue

        # Move the vertical line
        rightx = [newtimef, newtimef]
        righty = self.ui.mainplot.get_ylim()
        setp(self.rightslideline, xdata=rightx, ydata=righty)

        self.ui.graphicsView.draw()

        # Set the value to left slider
        self.ui.horizontalSlider_2.setValue(self._rightSlideValue)

        # Reset to line edit
        if resetT:
            self.ui.lineEdit_4.setText(str(newtimef))

        return

    def move_lowerSlider(self):
        """ Re-setup upper range line in figure.
        Triggered by a change in Qt Widget.  NO EVENT is required.
        """
        inewy = self.ui.verticalSlider_2.value()
        debug_msg = "LowerSlFider is set with value %s  vs. class variable %s" % (str(inewy), str(self._lowerSlideValue))
        Logger("Filter_Events").debug(debug_msg)

        # Return with no change
        if inewy == self._lowerSlideValue:
            # No change
            return

        if inewy >= self._upperSlideValue:
            # Out of upper range
            inewy = self._upperSlideValue - 1

        if inewy == 0 and self._lowerSlideValue < 0:
            setLineEdit = False
        else:
            setLineEdit = True

        # Move the lower vertical bar
        ylim = self.ui.mainplot.get_ylim()
        newy = ylim[0] + inewy*(ylim[1] - ylim[0])*0.01
        lowerx = self.ui.mainplot.get_xlim()
        lowery = [newy, newy]
        setp(self.lowerslideline, xdata=lowerx, ydata=lowery)

        self.ui.graphicsView.draw()

        # Set line edit input
        if setLineEdit is True:
            # Change value to line edit (5)
            self.ui.lineEdit_5.setText(str(newy))
            # Reset the class variable
            self._lowerSlideValue = inewy

        return

    def set_min_log_value(self):
        """ set the starting
        """
        debug_msg = "Minimum Log Value = %s" %(str(self.ui.lineEdit_minValue.text()))
        Logger("Filter_Events").debug(debug_msg)

        # get the limit of the canvas
        y_limits = self.ui.mainplot.get_ylim()

        if str(self.ui.lineEdit_minValue.text()) == "":
            # Empty. Default to minY
            new_y_min = y_limits[0]
        else:
            # Non empty.  Parse the input to float
            new_y_min = float(str(self.ui.lineEdit_minValue.text()))

        # Convert to integer slide value
        int_min_log_value = int( (new_y_min-y_limits[0])/(y_limits[1] - y_limits[0])*100 )
        debug_msg = "ilowerSlide = %s" % str(int_min_log_value)
        Logger("Filter_Events").debug(debug_msg)

        # Return if no change
        if int_min_log_value == self._lowerSlideValue:
            return

        # Set value if out of range
        resetL = True
        if int_min_log_value >= self._upperSlideValue:
            int_min_log_value = self._upperSlideValue - 1
        else:
            resetL = False

        if resetL is True:
            new_y_min = y_limits[0] + int_min_log_value * (y_limits[1]-y_limits[0]) * 0.01

        # Move the vertical line
        lowerx =  self.ui.mainplot.get_xlim()
        lowery =  [new_y_min, new_y_min]
        setp(self.lowerslideline, xdata=lowerx, ydata=lowery)

        self.ui.graphicsView.draw()

        # Move the slide bar (lower)
        self._lowerSlideValue = int_min_log_value
        debug_msg = "LineEdit5 set slide to %s" % str(self._lowerSlideValue)
        Logger("Filter_Events").debug(debug_msg)
        self.ui.verticalSlider_2.setValue(self._lowerSlideValue)

        # Reset line Edit if using default
        if resetL is True:
            self.ui.lineEdit_minValue.setText(str(new_y_min))

        return

    def move_upper_slider(self):
        """ Re-setup upper range line in figure.
        Triggered by a change in Qt Widget.  NO EVENT is required.
        """
        int_new_max_y = self.ui.verticalSlider.value()

        # Return w/o change
        if int_new_max_y == self._upperSlideValue:
            return

        # Set to boundary value
        if int_new_max_y <= self._lowerSlideValue:
            int_new_max_y = self._lowerSlideValue + 1

        # Reset line editor?
        if int_new_max_y == 100 and self._upperSlideValue > 100:
            to_set_y_max = False
        else:
            to_set_y_max = True

        # Move the upper value bar: upperx and uppery are real value (float but not (0,100)) of the figure
        ylim = self.ui.mainplot.get_ylim()
        newy = ylim[0] + int_new_max_y*(ylim[1] - ylim[0])*0.01
        upperx = self.ui.mainplot.get_xlim()
        uppery = [newy, newy]
        setp(self.upperslideline, xdata=upperx, ydata=uppery)

        self.ui.graphicsView.draw()

        # Change value to line editor max log value
        if to_set_y_max is True:
            self.ui.lineEdit_maxValue.setText(str(newy))
            self._upperSlideValue = int_new_max_y

        return

    def set_max_log_value(self):
        """ Set maximum log value from line-edit
        """
        # get value (in string) from line editor
        max_value_str = str(self.ui.lineEdit_maxValue.text())
        debug_msg = "Maximum Log Value = %s" % max_value_str
        Logger("Filter_Events").debug(debug_msg)

        canvas_y_limits = self.ui.mainplot.get_ylim()
        if len(max_value_str) == 0:
            # Empty. Default to minY
            new_max_y = canvas_y_limits[1]
        else:
            # Parse
            new_max_y = float(max_value_str)

        # synchronize to integer slide value
        slider_max_value = int((new_max_y-canvas_y_limits[0])/(canvas_y_limits[1] - canvas_y_limits[0])*100)
        debug_msg = "iUpperSlide = %s" % str(slider_max_value)
        Logger("Filter_Events").debug(debug_msg)

        if slider_max_value == self._upperSlideValue:
            # return if no change
            return
        else:
            # set to default if out of range
            resetL = True
            if slider_max_value < self._lowerSlideValue:
                slider_max_value = self._lowerSlideValue + 1
            else:
                resetL = False

            # Set newmaxY if necessary
            if resetL is True:
                new_max_y = canvas_y_limits[0] + slider_max_value * (canvas_y_limits[1] - canvas_y_limits[0]) * 0.01

            # Move the vertical line
            upperx =  self.ui.mainplot.get_xlim()
            uppery =  [new_max_y, new_max_y]
            setp(self.upperslideline, xdata=upperx, ydata=uppery)
            self.ui.graphicsView.draw()

            # Set the value to upper slider
            self._upperSlideValue = slider_max_value
            self.ui.verticalSlider.setValue(self._upperSlideValue)

            # set the value to editor if necessary
            if resetL:
                print '[DB...BAT] the value from the line editor is set back to line editor from some pamameters ' \
                      'from the canvas.'
                self.ui.lineEdit_maxValue.setText(str(new_max_y))
            # END-IF-resetL
        # END-IF

        return

    def browse_File(self):
        """ Open a file dialog to get file
        """
        filename = QtGui.QFileDialog.getOpenFileName(self, 'Input File Dialog',
                                                     self._defaultdir, "Data (*.nxs *.dat);;All files (*)")

        self.ui.lineEdit.setText(str(filename))

        info_msg = "Selected file: %s." % str(filename)
        Logger("Filter_Events").information(info_msg)

        return

    def load_file(self):
        """ Load the file by file name or run number
        """
        # Get file name from line editor
        file_name = str(self.ui.lineEdit_run.text())

        # get workspace
        data_ws = self._loadFile(str(file_name))
        if data_ws is None:
            error_msg = "Unable to locate run %s in default directory %s." % (file_name, self._defaultdir)
            Logger("Filter_Events").error(error_msg)
            self._setErrorMsg(error_msg)
        else:
            self._import_data_workspace(data_ws)
            self._defaultdir = os.path.dirname(str(file_name))

        # Reset GUI
        self._resetGUI(resetfilerun=False)

        return

    def use_existWS(self):
        """ Set up workspace to an existing one
        """
        wsname = str(self.ui.comboBox.currentText())

        try:
            dataws = AnalysisDataService.retrieve(wsname)
            self._import_data_workspace(dataws)
        except KeyError:
            pass

        # Reset GUI
        self._resetGUI(resetfilerun=True)

        return

    def plotLogValue(self):
        """ Plot log value
        """
        # Get log value
        log_name = str(self.ui.comboBox_sampleLogs.currentText())
        if len(log_name) == 0:
            # return due to the empty one is chozen
            return

        sample_property = self._dataWS.getRun().getProperty(log_name)
        vectimes = sample_property.times
        vecvalue = sample_property.value

        # check
        if len(vectimes) == 0:
            error_msg = "Empty log!"
            Logger("Filter_Events").error(error_msg)

        # convert absolute time to relative time in seconds
        t0 = self._dataWS.getRun().getProperty("proton_charge").times[0]
        t0ns = t0.totalNanoseconds()

        # append 1 more log if original log only has 1 value
        tf = self._dataWS.getRun().getProperty("proton_charge").times[-1]
        vectimes.append(tf)
        vecvalue = numpy.append(vecvalue, vecvalue[-1])

        vecreltimes = []
        for t in vectimes:
            rt = float(t.totalNanoseconds() - t0ns) * 1.0E-9
            vecreltimes.append(rt)

        # Set to plot
        xlim = [min(vecreltimes), max(vecreltimes)]
        ylim = [min(vecvalue), max(vecvalue)]
        self.ui.mainplot.set_xlim(xlim[0], xlim[1])
        self.ui.mainplot.set_ylim(ylim[0], ylim[1])

        setp(self.mainline, xdata=vecreltimes, ydata=vecvalue)

        samunit = sample_property.units
        if len(samunit) == 0:
            ylabel = log_name
        else:
            ylabel = "%s (%s)" % (log_name, samunit)
        self.ui.mainplot.set_ylabel(ylabel, fontsize=13)

        # assume that all logs are on almost same X-range.  Only Y need to be reset
        setp(self.leftslideline, ydata=ylim)
        setp(self.rightslideline, ydata=ylim)

        # reset the log value limit as previous one does not make any sense
        setp(self.lowerslideline, xdata=xlim, ydata=[ylim[0], ylim[0]])
        self._lowerSlideValue = 0
        self.ui.verticalSlider_2.setValue(self._lowerSlideValue)
        self.ui.lineEdit_minValue.setText("")

        setp(self.upperslideline, xdata=xlim, ydata=[ylim[1], ylim[1]])
        self._upperSlideValue = 100
        self.ui.verticalSlider.setValue(self._upperSlideValue)
        self.ui.lineEdit_maxValue.setText("")

        self.ui.graphicsView.draw()

        # Load property's statistic and give suggestion on parallel and fast log
        timeavg = sample_property.timeAverageValue()
        numentries = sample_property.size()
        stat = sample_property.getStatistics()

        duration = stat.duration
        mean = stat.mean
        freq = float(numentries)/float(duration)

        self.ui.label_mean.show()
        self.ui.label_meanvalue.show()
        self.ui.label_avg.show()
        self.ui.label_timeAvgValue.show()
        self.ui.label_freq.show()
        self.ui.label_freqValue.show()
        self.ui.label_logname.show()
        self.ui.label_lognamevalue.show()
        self.ui.label_logsize.show()
        self.ui.label_logsizevalue.show()

        self.ui.label_meanvalue.setText("%.5e"%(mean))
        self.ui.label_timeAvgValue.setText("%.5e"%(timeavg))
        self.ui.label_freqValue.setText("%.5e"%(freq))
        self.ui.label_lognamevalue.setText(log_name)
        self.ui.label_logsizevalue.setText(str(numentries))

        # Set suggested processing scheme
        if numentries > HUGE_FAST:
            self.ui.checkBox_fastLog.setCheckState(True)
            if numentries > HUGE_PARALLEL:
                self.ui.checkBox_doParallel.setCheckState(True)
            else:
                self.ui.checkBox_doParallel.setCheckState(False)
        else:
            self.ui.checkBox_fastLog.setCheckState(False)
            self.ui.checkBox_doParallel.setCheckState(False)

        return

    def _import_data_workspace(self, dataws):
        """ Import data workspace for filtering
        """
        if dataws is None:
            return

        # Plot time counts
        errmsg = self._plotTimeCounts(dataws)
        if errmsg is not None:
            errmsg = "Workspace %s has invalid sample logs for splitting. Loading \
                    failure! \n%s\n" % (str(dataws), errmsg)
            self._setErrorMsg(errmsg)
            return False

        # Import log
        self._sampleLogNames = [""]

        run = dataws.getRun()
        plist = run.getProperties()
        for p in plist:
            pv = p.value
            if isinstance(pv, numpy.ndarray):
                times = p.times
                if len(times) > 1:
                    self._sampleLogNames.append(p.name)
        # ENDFOR(p)

        # Set up sample log
        self.ui.comboBox_sampleLogs.clear()
        self.ui.comboBox_sampleLogs.addItems(self._sampleLogNames)

        # Side information
        self.ui.label_mean.hide()
        self.ui.label_meanvalue.hide()
        self.ui.label_avg.hide()
        self.ui.label_timeAvgValue.hide()
        self.ui.label_freq.hide()
        self.ui.label_freqValue.hide()

        # Hide 'log name' above the graphic view
        self.ui.label_logname.hide()
        self.ui.label_lognamevalue.hide()

        # Set dataws to class variable
        self._dataWS = dataws

        return True

    def scan_event_worksapces(self):
        """
        """
        wsnames = AnalysisDataService.getObjectNames()

        eventwsnames = []
        for wsname in wsnames:
            wksp = AnalysisDataService.retrieve(wsname)
            if wksp.__class__.__name__.count("Event") == 1:
                eventwsnames.append(wsname)
        # ENDFOR

        if len(eventwsnames) > 0:
            self.ui.comboBox.clear()
            self.ui.comboBox.addItems(eventwsnames)

        return

    def _loadFile(self, filename):
        """ Load file or run
        File will be loaded to a workspace shown in MantidPlot
        """
        config = ConfigService

        # Check input file name and output workspace name
        if filename.isdigit() is True:
            # Construct a file name from run number
            runnumber = int(filename)
            if runnumber <= 0:
                error_msg = "Run number cannot be less or equal to zero.  User gives %s. " % (filename)
                Logger("Filter_Events").error(error_msg)
                return None
            else:
                ishort = config.getInstrument(self._instrument).shortName()
                filename = "%s_%s" %(ishort, filename)
                wsname = filename + "_event"

        elif filename.count(".") > 0:
            # A proper file name
            wsname = os.path.splitext(os.path.split(filename)[1])[0]

        elif filename.count("_") == 1:
            # A short one as instrument_runnumber
            iname = filename.split("_")[0]
            str_runnumber = filename.split("_")[1]
            if str_runnumber.isdigit() is True and int(str_runnumber) > 0:
                # Acccepted format
                ishort = config.getInstrument(iname).shortName()
                wsname = "%s_%s_event" % (ishort, str_runnumber)
            else:
                # Non-supported
                error_msg = "File name / run number in such format %s is not supported. " % (filename)
                Logger("Filter_Events").error(error_msg)

                return None

        else:
            # Unsupported format
            error_msg = "File name / run number in such format %s is not supported. " % (filename)
            Logger("Filter_Events").error(error_msg)

            return None

        # Load
        try:
            ws = api.Load(Filename=filename, OutputWorkspace=wsname)
        except RuntimeError as e:
            ws = None
            return str(e)

        return ws

    def _plotTimeCounts(self, wksp):
        """ Plot time/counts
        """
        import datetime
        # Rebin events by pulse time
        try:
            # Get run start and run stop
            if wksp.getRun().hasProperty("run_start"):
                runstart = wksp.getRun().getProperty("run_start").value
            else:
                runstart = wksp.getRun().getProperty("proton_charge").times[0]
            runstop = wksp.getRun().getProperty("proton_charge").times[-1]

            runstart = str(runstart).split(".")[0].strip()
            runstop = str(runstop).split(".")[0].strip()

            t0 = datetime.datetime.strptime(runstart, "%Y-%m-%dT%H:%M:%S")
            tf = datetime.datetime.strptime(runstop, "%Y-%m-%dT%H:%M:%S")

            # Calcualte
            dt = tf-t0
            timeduration = dt.days*3600*24 + dt.seconds

            timeres = float(timeduration)/MAXTIMEBINSIZE
            if timeres < 1.0:
                timeres = 1.0

            sumwsname = "_Summed_%s"%(str(wksp))
            if AnalysisDataService.doesExist(sumwsname) is False:
                sumws = api.SumSpectra(InputWorkspace=wksp, OutputWorkspace=sumwsname)
                sumws = api.RebinByPulseTimes(InputWorkspace=sumws, OutputWorkspace = sumwsname,
                                              Params="%f"%(timeres))
                sumws = api.ConvertToPointData(InputWorkspace=sumws, OutputWorkspace=sumwsname)
            else:
                sumws = AnalysisDataService.retrieve(sumwsname)
        except RuntimeError as e:
            return str(e)

        vecx = sumws.readX(0)
        vecy = sumws.readY(0)

        xmin = min(vecx)
        xmax = max(vecx)
        ymin = min(vecy)
        ymax = max(vecy)

        # Reset graph
        self.ui.mainplot.set_xlim(xmin, xmax)
        self.ui.mainplot.set_ylim(ymin, ymax)

        self.ui.mainplot.set_xlabel('Time (seconds)', fontsize=13)
        self.ui.mainplot.set_ylabel('Counts', fontsize=13)

        # Set up main line
        setp(self.mainline, xdata=vecx, ydata=vecy)

        # Reset slide
        newslidery = [min(vecy), max(vecy)]

        newleftx = xmin + (xmax-xmin)*self._leftSlideValue*0.01
        setp(self.leftslideline, xdata=[newleftx, newleftx], ydata=newslidery)

        newrightx = xmin + (xmax-xmin)*self._rightSlideValue*0.01
        setp(self.rightslideline, xdata=[newrightx, newrightx], ydata=newslidery)

        self.ui.graphicsView.draw()

        return

    def filter_by_time(self):
        """ Filter by time
        """
        # Generate event filters
        kwargs = {}
        if str(self.ui.lineEdit_startTime.text()) != "":
            rel_starttime = float(self.ui.lineEdit_startTime.text())
            kwargs["StartTime"] = str(rel_starttime)
        if self.ui.lineEdit_stopTime.text() != "":
            rel_stoptime = float(self.ui.lineEdit_stopTime.text())
            kwargs["StopTime"] = str(rel_stoptime)
        if self.ui.lineEdit_timeInterval.text() != "":
            interval = float(self.ui.lineEdit_timeInterval.text())
            kwargs["TimeInterval"] = interval

        splitwsname = str(self._dataWS) + "_splitters"
        splitinfowsname = str(self._dataWS) + "_info"

        title = str(self.ui.lineEdit_title.text())
        fastLog = self.ui.checkBox_fastLog.isChecked()

        splitws, infows = api.GenerateEventsFilter(
            InputWorkspace      = self._dataWS,
            UnitOfTime          = "Seconds",
            TitleOfSplitters    = title,
            OutputWorkspace     = splitwsname,
            FastLog             = fastLog,
            InformationWorkspace = splitinfowsname, **kwargs)

        self.splitWksp(splitws, infows)

        return

    def filter_by_log_value(self):
        """ Filter by log value
        """
        # Generate event filter
        kwargs = {}
        samplelog = str(self.ui.comboBox_sampleLogs.currentText())
        if len(samplelog) == 0:
            error_msg = "No sample log is selected!"
            Logger("Filter_Events").error(error_msg)
            return

        if str(self.ui.lineEdit_startTime.text()) != '':
            rel_starttime = float(self.ui.lineEdit_startTime.text())
            kwargs["StartTime"] = str(rel_starttime)

        if self.ui.lineEdit_stopTime.text() != "":
            rel_stoptime = float(self.ui.lineEdit_stopTime.text())
            kwargs["StopTime"] = str(rel_stoptime)

        if self.ui.lineEdit_minValue.text() != "":
            minlogvalue = float(self.ui.lineEdit_minValue.text())
            kwargs["MinimumLogValue"] = minlogvalue

        if self.ui.lineEdit_maxValue.text() != "":
            maxlogvalue = float(self.ui.lineEdit_maxValue.text())
            kwargs["MaximumLogValue"] = maxlogvalue

        if self.ui.lineEdit_logStep.text() != "":
            logvalueintv = float(self.ui.lineEdit_logStep.text())
            kwargs["LogValueInterval"] = logvalueintv
        logvalchangedir = str(self.ui.comboBox_4.currentText())
        kwargs["FilterLogValueByChangingDirection"] = logvalchangedir

        if self.ui.lineEdit_9.text() != "":
            logvalueintv = float(self.ui.lineEdit_9.text())
            kwargs["TimeTolerance"] = logvalueintv
        logboundtype = str(self.ui.comboBox_5.currentText())
        kwargs["LogBoundary"] = logboundtype

        if self.ui.lineEdit_8.text() != "":
            logvaluetol = float(self.ui.lineEdit_8.text())
            kwargs["LogValueTolerance"] = logvaluetol

        splitwsname = str(self._dataWS) + "_splitters"
        splitinfowsname = str(self._dataWS) + "_info"
        fastLog = self.ui.checkBox_fastLog.isChecked()

        title = str(self.ui.lineEdit_title.text())

        splitws, infows = api.GenerateEventsFilter(
            InputWorkspace      = self._dataWS,
            UnitOfTime          = "Seconds",
            TitleOfSplitters    = title,
            OutputWorkspace     = splitwsname,
            LogName             = samplelog,
            FastLog             = fastLog,
            InformationWorkspace = splitinfowsname, **kwargs)

        try:
            self.splitWksp(splitws, infows)
        except RuntimeError as e:
            self._setErrorMsg("Splitting Failed!\n %s" % (str(e)))

        return

    def splitWksp(self, splitws, infows):
        """ Run FilterEvents
        """
        dogroupws = self.ui.checkBox_groupWS.isChecked()
        filterbypulse = self.ui.checkBox_filterByPulse.isChecked()
        startfrom1 = self.ui.checkBox_from1.isChecked()
        splitsamplelog = self.ui.checkBox_splitLog.isChecked()

        corr2sample = str(self.ui.comboBox_tofCorr.currentText())
        how2skip = str(self.ui.comboBox_skipSpectrum.currentText())

        kwargs = {}
        if corr2sample == "Direct":
            ei = float(self.ui.lineEdit_Ei.text())
            kwargs["IncidentEnergy"] = ei
        elif corr2sample == "Customized":
            corrws = str(self.ui.comboBox_corrWS.currentText())
            kwargs["DetectorTOFCorrectionWorkspace"] = corrws

        # Output workspace name
        outbasewsname = str(self.ui.lineEdit_outwsname.text())
        if len(outbasewsname) == 0:
            outbasewsname = "tempsplitted"
            self.ui.lineEdit_outwsname.setText(outbasewsname)

        api.FilterEvents(
            InputWorkspace          = self._dataWS,
            SplitterWorkspace       = splitws,
            InformationWorkspace    = infows,
            OutputWorkspaceBaseName = outbasewsname,
            GroupWorkspaces         = dogroupws,
            FilterByPulseTime       = filterbypulse,
            CorrectionToSample      = corr2sample,
            SpectrumWithoutDetector = how2skip,
            SplitSampleLogs         = splitsamplelog,
            OutputWorkspaceIndexedFrom1     = startfrom1,
            OutputTOFCorrectionWorkspace    = 'TOFCorrTable', **kwargs)

        return

    def showHideEi(self):
        """
        """
        corrtype = str(self.ui.comboBox_tofCorr.currentText())

        # Incident energy
        if corrtype == "Direct":
            self.ui.label_Ei.show()
            self.ui.lineEdit_Ei.show()
        else:
            self.ui.label_Ei.hide()
            self.ui.lineEdit_Ei.hide()

        # Workspace
        if corrtype == "Customized":
            self.ui.label_Ei_2.show()
            self.ui.comboBox_corrWS.show()
            self.ui.pushButton_refreshCorrWSList.show()

            # Search for table workspace
            self._searchTableWorkspaces()

        else:
            self.ui.label_Ei_2.hide()
            self.ui.comboBox_corrWS.hide()
            self.ui.pushButton_refreshCorrWSList.hide()

        return

    def _searchTableWorkspaces(self):
        """ Search table workspaces and add to 'comboBox_corrWS'
        """
        wsnames = AnalysisDataService.getObjectNames()

        tablewsnames = []
        for wsname in wsnames:
            wksp = AnalysisDataService.retrieve(wsname)
            if isinstance(wksp, mantid.api.ITableWorkspace):
                tablewsnames.append(wsname)
        # ENDFOR

        self.ui.comboBox_corrWS.clear()
        if len(tablewsnames) > 0:
            self.ui.comboBox_corrWS.addItems(tablewsnames)

        return

    def _clearErrorMsg(self):
        """ Clear error message
        """
        #self.ui.plainTextEdit_ErrorMsg.setPlainText("")
        #self.ui.label_error.hide()

        return

    def _setErrorMsg(self, errmsg):
        """ Clear error message
        """
        #self.ui.plainTextEdit_ErrorMsg.setPlainText(errmsg)
        #self.ui.label_error.show()

        #print "Testing Pop-up Error Message Window: %s" % (errmsg)
        self._errMsgWindow = MyPopErrorMsg()
        self._errMsgWindow.setMessage(errmsg)
        self._errMsgWindow.show()

        return

    def helpClicked(self):
        from pymantidplot.proxies import showCustomInterfaceHelp
        showCustomInterfaceHelp("FilterEventUI")

    def _resetGUI(self, resetfilerun=False):
        """ Reset GUI including all text edits and etc.
        """
        if resetfilerun is True:
            self.ui.lineEdit.clear()

        # Plot related
        self.ui.lineEdit_3.clear()
        self.ui.lineEdit_4.clear()
        self.ui.horizontalSlider.setValue(0)
        self.ui.horizontalSlider_2.setValue(100)

        self.ui.lineEdit_outwsname.clear()
        self.ui.lineEdit_title.clear()

        # Filter by log value
        self.ui.lineEdit_5.clear()
        self.ui.lineEdit_6.clear()

        self.ui.verticalSlider_2.setValue(0)
        self.ui.verticalSlider.setValue(100)

        ylim = self.ui.mainplot.get_ylim()
        miny = ylim[0]
        maxy = ylim[1]
        xlim = self.ui.mainplot.get_xlim()
        setp(self.lowerslideline, xdata=xlim, ydata=[miny, miny])
        setp(self.upperslideline, xdata=xlim, ydata=[maxy, maxy])
        self.ui.graphicsView.draw()

        self.ui.lineEdit_7.clear()
        self.ui.lineEdit_8.clear()
        self.ui.lineEdit_9.clear()

        # Filter by time
        self.ui.lineEdit_timeInterval.clear()

        # Advanced setup
        self.ui.comboBox_tofCorr.setCurrentIndex(0)
        self.ui.lineEdit_Ei.clear()

        self.ui.checkBox_fastLog.setCheckState(False)
        self.ui.checkBox_doParallel.setCheckState(False)

        self.ui.comboBox_skipSpectrum.setCurrentIndex(0)

        self.ui.checkBox_filterByPulse.setCheckState(False)
        self.ui.checkBox_from1.setCheckState(False)
        self.ui.checkBox_groupWS.setCheckState(True)
        self.ui.checkBox_splitLog.setCheckState(False)

        # Error message
        # self.ui.plainTextEdit_ErrorMsg.clear()

        return
