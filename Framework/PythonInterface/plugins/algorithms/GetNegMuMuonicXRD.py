from __future__ import (absolute_import, division, print_function)
from mantid.api import * # PythonAlgorithm, registerAlgorithm, WorkspaceProperty
from mantid.kernel import *
from mantid.simpleapi import *

#pylint: disable=no-init


class GetNegMuMuonicXRD(PythonAlgorithm):
    #Dictionary of <element>:<peaks> easily extendible by user.
    muonic_xr ={'He' :[10.484,10.278,9.742,8.228], #[data point 2: 10.34,9.70,8.18] [data point 3:8.18]  
	            'Li' :[18.69],#[li'7 ,18.1]  [li'6: 18.64,18.1]  
				'Be' :[33.39],#[data point 2:33.0]
				'B'  :[52.31], #(nutron number 11[data point 2:52.23] [data point 3:51.6])
			     #'B' (Z = 10 [data point 1: 51.6] [data point 2:52.23] [data point 3:52.18])
				'C'  :[97.601,96.355,94.095,89.212,75.248],# [data point 2:75.8] [data point 3: 94.0, 88.6,75.3]
                 # [data point 4: 75.25] [data point 5:75.23]	
     			'N'  :[131.167,96.355,94.095,89.212,75.248],# [data point 2:101.9] [data point 3:102.29]
		    	'O'  :[173.331,171.144,167.114,158.408,133.525], #[data point 2:133.4] [data point 3:158.9,133.3] [data point 4:133.56]
				'F'  :[268.45], #[data point 2:168.9] [data point 3:168.07]
				'Ne' :[207.1],
				'Na' :[250.21], #[data point 2:249.6] [data point 3: 250.24]
				'Mg' :[382.00,372.50,353.00,296.40], # [data point 2: 296.1] [data point 3: 296.88] [data point 4:296.55]
				'Ar' :[452.62,446.66,436.14,413.02,346.90], #[data point 2:347.21] [data point 3:346.82]
				'Si' :[503.85,477.00,400.15], # [data point 2:400.38] [data point 3:400.22] [data point 4:400.2]
				'P'  :[456.54], # [data point 2:457.06]
				'S'  :[651.90,617.00,516.50], # [data point 2: 516.24] 
				'Cl' :[578.56],
				'Al' :[770.6,643.94],# (Z=36 [data point 1:644.34]) (Z=38 [644.34]) (z=40 [643.94]) 
				'K'  :[712.24], #[data point 1:712.64]
				'Ca' :[940.7,783.85,212.05,158.17,156.830],
				'Ti' :[1011.3,191.921,189.967],
				'V'  :[1094.4,208.1],
				'Cr' :[1094.4],
				'Mn' :[1171.2],
		 	    'Fe' :[1257.21,1253.01,269.39,265.69,264.95,244.33,243.17],
				'Co' :[1541.8],
				'Ni' :[1427.4,1422.1,309.97],
		    	'As' :[1866.9,1855.8,436.6,427.5],
	            'Au' :[8135.2,8090.6,8105.4,8069.4,5764.89,5594.97,3360.2,
                        3206.8,2474.22,2341.21,2304.44,1436.05,1391.58,1104.9,
                        899.14,869.98,405.654,400.143],
                'Ag' :[3184.7,3147.7,901.2,869.2,308.428,304.759],
                'Cu' :[1512.78,1506.61,334.8,330.26],
			    'Sr' :[2241.5,2224.0,200.254,198.708],
				'Y'  :[3038.63,3033.11,2439.38,2420.12,807.6,616.38,599.39,200.254,
				       198.708],
				'zr' :[2535.9,2514.9],
				'Nb' :[2626.6,2603.42,683.03,662.54,543.91,537.45,140.01,116.87],
				'Mo' :[2732.3,2706.8,717.8,695.0],
				'Rh' :[2982],
				'Pd' :[3077],
				'Cd' :[3263.63,3223.74,940.58,905.85,321.973,317.977],
				'In' :[3366.27,3322.67,981.64,943.39],
				'Sn' :[3454.41,3408.79,1369.50,1325.9,1022.2,982.20,976.55,747.15,
				       733.20,280.14,234.50,350.0,345.3],
				'Sb' :[3543.3,3497.7,1062.8,1019.6,360.3],
				'Te' :[3625.6,3375.5,1104.3,1060.0,375.9],
				'I'  :[3723.23,3667.35,1150.42,1101.84,394.30,388.16,179.98],
				'Cs' :[3899.1,3836.1,421.2],
				'Ba' :[5196.7,5215.3,3988.09,3922.15,1720.40,1658.31,1283.22,
				       1227.09,1218.03,868.50,405.41,339.5,441.28,433.72,
					   200.53],
				'La' :[4081.42,4011.95,1329.90,1270.13,1259.02,452.8,206.8],
				'Ce' :[5449.2,5469.6,5459.7,5459.7,4155.86,4082.81,1846.10,
				       1777.27,1375.80,1313.47,1303.10,932.55,913.00,453.45,
					   379.90,474.24,465.46,215.01],
				'Pr' :[4263.54,4185.87,1424.07,1158.11,134.43,485.5,221.3],
				'Nd' :[4352.52,4270.78,1472.30,1403.06,1390.52,498.4,229.4],
				'Hg' :[5817.2,5645.1,2523.6,2388.5,918.1,888.1,412.7],
				'Tl' :[5906.38,5726.14,2587.73,2448.22,2407.48,1485.67,
				       1439.54,1178.52,947.47,915.1,426.828,420.717],        
                 'Zn' :[1600.15,1592.97,360.75,354.29],
                 'Pb' :[8523.3,8442.11,5966.0,5780.1,2641.8,2499.7,
                        2459.7,1511.63,1214.12,1028.83,972.3,938.4,
                        437.687,431.285],
                'As' :[1866.9,1855.8,436.6,427.5],
                'Sn' :[3457.3,3412.8,1022.6,982.5,349.953,345.226],
				'Bi' :[6032.2,5839.7,2700.2,2554.8,2554.8,2501.81,996.53,
				       961.00,448.80,442.11],
				}

    def PyInit(self):
        self.declareProperty(StringArrayProperty("Elements", values=[],
                                                 direction=Direction.Input))
        self.declareProperty(name="YAxisPosition",
                             defaultValue=-0.001,
                             doc="Position for Markers on the y-axis")
        self.declareProperty(WorkspaceGroupProperty('OutputWorkspace', '', direction=Direction.Output),
                             doc='The output workspace will always be a GroupWorkspaces '
                                 'that will have the workspaces of each'
                                 ' muonicXR workspace created')
        #We sort the lists of x-values from the dictionary here
        #so that mantid can plot the workspaces it produces.
        for key in self.muonic_xr:
            value = self.muonic_xr.get(key)
            self.muonic_xr[key] = sorted(value)

    def get_muonic_xr(self, element):
        #retrieve peak values from dictionary Muonic_XR
        peak_values = self.muonic_xr[element]
        return peak_values

    def validateInput(self):
        issues = dict()

        elements = self.getProperty('Elements').value()
        if elements == "":
            issues["Elements"] = 'No elements have been selected from the periodic table'
        y_axis_position = self.getProperty('YAxisPosition').value()
        if y_axis_position == "":
            issues["YAxisPosition"] = 'No y-shift value has been entered'
        outputworkspace_str = self.getProperty('OutputWorkspace').value()
        if outputworkspace_str == "":
            issues['OutputWorkspace'] = 'No output workspace name has been specified'

        return issues

    def create_muonic_xr_ws(self, element, y_pos):
        #retrieve the values from Muonic_XR
        xr_peak_values = self.get_muonic_xr(element)
        #Calibrate y-axis for workspace
        y_pos_ws = [y_pos]*len(xr_peak_values)
        xvalues = xr_peak_values
        muon_xr_ws = CreateWorkspace(xvalues, y_pos_ws[:])
        RenameWorkspaces(muon_xr_ws, WorkspaceNames="MuonXRWorkspace_"+element)
        return muon_xr_ws

    def category(self):
        return "Muon"

    def PyExec(self):
        elements = self.getProperty("Elements").value
        y_position = self.getProperty("YAxisPosition").value
        workspace_list = [None]*len(elements)
        for idx,element in enumerate(elements):
            curr_workspace = self.create_muonic_xr_ws(element, y_position)
            workspace_list[idx] = curr_workspace

        self.setProperty("OutputWorkspace", GroupWorkspaces(workspace_list))
        self.log().information(str("Created Group: "+ self.getPropertyValue("OutputWorkspace")))

AlgorithmFactory.subscribe(GetNegMuMuonicXRD)
