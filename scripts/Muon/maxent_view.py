from __future__ import (absolute_import, division, print_function)

from PyQt4 import QtCore, QtGui

from Muon import table_utils


class MaxEntView(QtGui.QWidget):
    """
    The view for the MaxEnt widget. This
    creates the look of the widget
    """
    # signals
    maxEntButtonSignal = QtCore.pyqtSignal()
    cancelSignal = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        super(MaxEntView, self).__init__(parent)
        self.grid = QtGui.QVBoxLayout(self)

        #make table
        self.table = QtGui.QTableWidget(self)
        self.table.resize(800, 800)

        self.table.setRowCount(9)
        self.table.setColumnCount(2)
        self.table.setColumnWidth(0,300)
        self.table.setColumnWidth(1,300)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setHorizontalHeaderLabels(("MaxEnt Property;Value").split(";"))
        table_utils.setTableHeaders(self.table)

        # populate table
        options = ['test']

        table_utils.setRowName(self.table,0,"Workspace")
        self.ws = table_utils.addComboToTable(self.table,0,options)
 
        table_utils.setRowName(self.table,1,"First good time")
        self.first_good= table_utils.addDoubleToTable(self.table,0.1,1)
 
        table_utils.setRowName(self.table,2,"Last good time")
        self.last_good= table_utils.addDoubleToTable(self.table,15.0,2)

        table_utils.setRowName(self.table, 3,"Fit dead times")
        self.dead_box= table_utils.addCheckBoxToTable(self.table,True,3)

        table_utils.setRowName(self.table,4,"Construct Phase Table")
        self.phaseTable_box = table_utils.addCheckBoxToTable(self.table,True,4)

        table_utils.setRowName(self.table, 5,"Fix phases")
        self.fix_phase_box= table_utils.addCheckBoxToTable(self.table,False,5)
 
        table_utils.setRowName(self.table, 6,"Output phase table")
        self.output_phase_box= table_utils.addCheckBoxToTable(self.table,False,6)

        table_utils.setRowName(self.table, 7,"Output deadtimes")
        self.output_dead_box= table_utils.addCheckBoxToTable(self.table,False,7)

        table_utils.setRowName(self.table, 8,"Output reconstructed data")
        self.output_data_box= table_utils.addCheckBoxToTable(self.table,False,8)

 



        self.table.resizeRowsToContents()

        # advanced options table
        self.advancedLabel=QtGui.QLabel("\n  Advanced Options")
        #make table
        self.tableA = QtGui.QTableWidget(self)
        self.tableA.resize(800, 800)

        self.tableA.setRowCount(6)
        self.tableA.setColumnCount(2)
        self.tableA.setColumnWidth(0,300)
        self.tableA.setColumnWidth(1,300)

        self.tableA.verticalHeader().setVisible(False)
        self.tableA.horizontalHeader().setStretchLastSection(True)

        self.tableA.setHorizontalHeaderLabels(("Advanced Property;Value").split(";"))
        table_utils.setTableHeaders(self.tableA)

        table_utils.setRowName(self.tableA,0,"Maximum entropy constant (A)")
        self.AConst= table_utils.addDoubleToTable(self.tableA,0.1,0)
 
        table_utils.setRowName(self.tableA,1,"Lagrange multiplier for chi^2")
        self.factor= table_utils.addDoubleToTable(self.tableA,1.04,1)

        table_utils.setRowName(self.tableA,2,"Inner Iterations")
        self.inner_loop= table_utils.addSpinBoxToTable(self.tableA,10,2)

        table_utils.setRowName(self.tableA,3,"Outer Iterations")
        self.outerr_loop= table_utils.addSpinBoxToTable(self.tableA,10,3)

        table_utils.setRowName(self.tableA, 4,"Double pulse data")
        self.double_pules_box= table_utils.addCheckBoxToTable(self.tableA,False,4)

        table_utils.setRowName(self.tableA,5,"Number of data points")
        self.N_points = table_utils.addComboToTable(self.tableA,5,options)



        #layout
        # this is if complex data is unhidden
        self.table.setMinimumSize(40,203)
        self.tableA.setMinimumSize(40,207)
        self.horizontalSpacer1 = QtGui.QSpacerItem(20, 30, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        self.horizontalSpacer2 = QtGui.QSpacerItem(20, 70, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        #make buttons
        self.button = QtGui.QPushButton('Calculate MaxEnt', self)
        self.button.setStyleSheet("background-color:lightgrey")
        self.cancel = QtGui.QPushButton('Cancel', self)
        self.cancel.setStyleSheet("background-color:lightgrey")
        self.cancel.setEnabled(False)
        #connects
        self.button.clicked.connect(self.MaxEntButtonClick)
        self.cancel.clicked.connect(self.cancelClick)
        # button layout
        self.buttonLayout=QtGui.QHBoxLayout()
        self.buttonLayout.addWidget(self.button)
        self.buttonLayout.addWidget(self.cancel)
        # add to layout
        self.grid.addWidget(self.table)
        self.grid.addItem(self.horizontalSpacer1)
        self.grid.addWidget(self.advancedLabel)
        self.grid.addWidget(self.tableA)
        self.grid.addItem(self.horizontalSpacer2)
        self.grid.addLayout(self.buttonLayout)

    # add data to view
    def addItems(self,options):
        self.ws.clear()
        self.ws.addItems(options)

    # send signal
    def MaxEntButtonClick(self):
        self.maxEntButtonSignal.emit()

    def cancelClick(self):
        self.cancelSignal.emit()

    # get some inputs for model
    def initMaxEntInput(self):
        inputs={}

        #  this will be removed once maxEnt does a simultaneous fit
        inputs['InputWorkspace']=str( self.ws.currentText()).replace(";","; ")
        # will use this instead of the above
        #inputs['InputWorkspace']="MuonAnalysis"
        inputs['ComplexData']=  self.complex_data_box.checkState()
        inputs["ComplexImage"]=  self.complex_image_box.checkState()
        inputs['PositiveImage']=self.positive_image_box.checkState()
        inputs["ResolutionFactor"]=int(self.resolution_box.text())
        inputs["A"] = float(self.AConst.text())
        inputs["AutoShift"]=self.shift_box.checkState()
        inputs["ChiTargetOverN"]=float(self.chiTarget.text())
        inputs["ChiEps"]=float(self.chiEps.text())
        inputs["DistancePenalty"]=float(self.dist.text())
        inputs["MaxAngle"]=float(self.angle.text())
        inputs["MaxIterations"]=int(self.max_iterations.text())
        inputs["AlphaChopIterations"]=int(self.chop.text())

        # will remove this when sim maxent Works
        out=str( self.ws.currentText()).replace(";","; ")

        inputs['EvolChi']=out+";EvolChi;MaxEnt"
        inputs['EvolAngle']=out+";EvolAngle;MaxEnt"
        inputs['ReconstructedImage']=out+";FrequencyDomain;MaxEnt"
        inputs['ReconstructedData']=out+";TimeDomain;MaxEnt"

        return inputs

    def addRaw(self,inputs,key):
        inputs[key]+="_Raw"

    def isRaw(self):
        return self.raw_box.checkState() == QtCore.Qt.Checked

    # turn button on and off
    def activateCalculateButton(self):
        self.button.setEnabled(True)
        self.cancel.setEnabled(False)

    def deactivateCalculateButton(self):
        self.button.setEnabled(False)
        self.cancel.setEnabled(True)
