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
				

				#below is off a graph so not as acurate 
				'C'  :[60,100],
				'N'  :[101.9,121,128.1,131.2],
				'O'  :[18,26,31,76,135,158,178],
				'F'  :[109,175,215,230,267,539,915,1623,2051],
				'Ne' :[95,215,305,545],
				'Na' :[190,280,320,1285,1970],
				'Mg' :[100,310,460],
				'Al' :[90,400,480,1115,1860],				
				'Si' :[830,1025,1830],	
				'P'  :[130,580,2240],		
				'S'  :[510,1305,2260],		
				'Cl' :[140,480,700,785,1285,2130],
				'Ar' :[93,102,130,172,207,340],
				'K'  :[120,205,740,880,955,1615,2170],
				'Ca' :[1075,1285,2190,2515,2720,3010],
				'Ti' :[100,275,925,1150,1280],
				'V'  :[1050,1230,1350,1595],
				'Cr' :[120,320,1085,1265,1475],
				'Mn' :[70,305,1405,1610],
				'Fe' :[120,280,1265,1565,1755],
				'Co' :[270,1730,1860],
				'Ni' :[120,190,315,400,860,1445,1730,1990],
				'As' :[100,425,605,1790,2300,2620],
				'Au' :[130,380,595,890,2360,2490,4740,5080,5270,5670,5860],
				'Ag' :[85,170,340,530,870,1150,2185,3210,4040,4750],
				'Cu' :[130,335,410,1140,1630,1820,2100],
				'Sr' :[150,315,540,760,1850,2360,2920,3345],
				'Y'  :[110,245,585,865,1940,2445,3035,3455],
				'Zr' :[95,180,340,660,860,1500,1730,2075,2445,3150,3640],
				'Nb' :[110,240,340,630,960,1645,2160,2655,3320,3860,5130],
				'Mo' :[135,260,345,785,930,1645,2145,2740,3460,3985],
				'Rh' :[95,130,295,800,1110,1205,1430,1980,2485,2990,3250,3640,3900,4410],
				'Pd' :[120,285,840,1175,2135,2555,3115,3930,4565],
				'Cd' :[165,235,330,505,920,1265,1355,2200,2650,3280,4190,4850],
				'In' :[90,170,365,645,935,1265,2375,3285,4315,5050],
				'Sn' :[160,365,520,1000,1350,2470,2980,3465,4560,5240],
				'Sb' :[110,210,375,530,1055,1430,2450,3030,3540,4040,4640],
				'Te' :[130,195,380,1090,1450,2510,3065,3680,4605],
				'I'  :[145,185,255,405,690,1160,1450,2785,3230,3800,4835],
				'Cs' :[90,230,440,585,1195,1650,2930,3350,3920,4575,5095,5510],
				'Ba' :[105,450,1250,1670,2910,3480,3965,4630,5210,6260],
				'La' :[90,195,475,1310,1400,1740,3030,3500,3960,4850,5375],
				'Ce' :[100,300,520,730,1430,1845,3030,3620,4120,4945,5500,6000],
				'Pr' :[135,245,490,560,1355,1640,3185,4165,5130,5640],
				'Nd' :[170,245,515,1405,1500,1990,3740,4245,4350,5745],
				'Hg' :[135,235,985,2400,2560,4640,5250,5425,5645,5870],
				'Tl' :[145,255,460,935,1330,2045,2450,2600,3400,4710,5180,5900],
				'Zn' :[50,110,330,460,1640,1950,2235],
				'Pb' :[90,420,980,1350,1940,2540,4975,5350,5530,5970],
				'Bi' :[100,240,495,975,2570,2710,3435,4670,5100,5325,5605,6045],
							
				'Sc' :[78,129,154,560,890,1010,1175,1250,1355,1525,2160],
				'Ga' :[105,180,320,480,530,600,910,1615,1790,1870,2090,2200,2230,2350],
                'Ge' :[85,100,130,190,210,410,500,540,660,700,1775,2155,2295,2345,2460],
                'Se' :[95,130,180,225,250,470,545,610,690,1440,1985,2320,2400,2490,2655,2715],
				'Br' :[445,540,615,730,1035,1440,2020,2210,2620,2705,2810,2900],
				'Rb' :[105,160,260,305,535,725,815,865,930,1235,1715,2340,2770,2965,3190],
				'Ru' :[90,125,295,760,1015,1140,1350,1880,2345,2905,3150,2660,3130,3760,3260],
				'Cs' :[100,200,295,355,580,620,690,730,1165,1230,1600,1655,1795,2810,2900,3330,3400,3845,3905],
				'Ba' :[75,140,205,300,400,505,635,1210,1290,1645,1740,1805,1880,2895
				       ,2980,3405,3490,3920,3995,4775,5215,5645,6195],
				'La' :[450,500,575,770,124,1310,400,1705,1790,1900,2230,2630,2980,3055
				       ,3475,3565,4010,4860,4840,5350,5830,6380],
				'Hf' :[90,125,260,320,520,715,1050,1190,1980,2120,2640,2890,4030,4235,4885,5020,5240],
				'Ta' :[90,190,300,345,505,745,1080,1230,1905,1955,2015,2030,2185,2950,4110,4720,5140,5315],
				'W'  :[115,235,330,560,600,765,1100,2025,2810,2990,4105,4610,5120,5230
				       ,5395,6550,6960,7450,8650,9145],
				'Re' :[105,160,220,325,420,800,1140,1740,2095,2180,3070,4260,4885,5295,5460,6555,7050,7630],
				'Os' :[60,130,285,815,1000,1160,2115,2300,4350,4850,4985,5460,5600],
				'Ir' :[110,140,210,390,600,890,1210,1400,1885,2120,2280,2370,2980,3120,4405,4595,4945,5075,7310,7850],
				'Pt' :[90,220,315,380,870,1250,1445,1800,2300,2400,2650,3260,3300,4600,4675,4965,5000,5160,5500,5700,7450,7965],
				'Sm' :[110,200,320,360,555,800,1545,1605,2010,3390,3475,3860,3980,4390,4510,5450,6000,6505,6635,7175],
				'Eu' :[100,160,305,535,730,1000,1480,1560,1650,2100,2740,3420,3900,4400,4465,4685,5550,6090],
				'Gd' :[90,190,430,530,830,1550,3990,4040,4445,4500,4610,5680,6195],
				'Tb' :[65,100,150,295,500,865,1590,1700,2305,3470,4075,4565,4590,4740,5330,5810,6300,7900],
				'Ho' :[100,175,290,595,640,745,900,1670,1870,2485,3695,4290,4565,4750,4900],
				'Er' :[150,290,525,965,1860,1990,3450,4290,4870,4930,6205,6745],
				'Tm' :[105,185,480,760,970,1800,1905,2435,3860,4000,4375,4830,5000,6470,6820],
				'Yb' :[115,130,365,505,755,1000,1865,2000,2555,3970,4120,4460,4895,5155,6440],
				'Lu' :[105,245,315,715,1020,1840,2000,2075,2770,3990,4555,5000,5050,5145,6450,7045],
				'Th' :[125,145,335,405,810,1100,1200,17352,900,3215,4300,5030,5530,6090,6300,6350],

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
