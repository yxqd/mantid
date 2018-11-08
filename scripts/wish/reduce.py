import numpy as n
import sys

import mantid.simpleapi as mantid
from mantid.simpleapi import *

sys.path.insert(0, "/opt/Mantid/bin")
sys.path.append("/isis/NDXWISH/user/scripts/autoreduction")


class Wish:

    def __init__(self, input_mode, cal_directory, user_directory, output_folder, delete_workspace):
        self.name = input_mode
        self.cal_dir = cal_directory
        self.use_folder = user_directory
        self.out_folder = output_folder
        self.deleteWorkspace = delete_workspace
        self.username = None
        self.data_directory = None
        self.user_directory = None
        self.datafile = None
        self.userdatadir = None
        self.userdataprocessed = None
        self.return_panel = {
            1: (6, 19461),
            2: (19462, 38917),
            3: (38918, 58373),
            4: (58374, 77829),
            5: (77830, 97285),
            6: (97286, 116741),
            7: (116742, 136197),
            8: (136198, 155653),
            9: (155654, 175109),
            10: (175110, 194565),
            0: (6, 194565)
        }

    def validate(self):
        """
        Autoreduction validate Function
        -------------------------------

        Function to ensure that the files we want to use in reduction exist.
        Please add any files/directories to the required_files/dirs lists.
        """
        print("Running validation")
        required_files = [self.datafile]
        required_dirs = [self.user_directory]
        for file_path in required_files:
            if not os.path.isfile(file_path):
                raise RuntimeError("Unable to find file: {}".format(file_path))
        for folder in required_dirs:
            if not os.path.isdir(folder):
                raise RuntimeError("Unable to find directory: {}".format(folder))
        print("Validation successful")

    def set_user_name(self, username):
        self.username = username

    def set_data_directory(self, directory="/archive/ndxwish/Instrument/data/cycle_09_5/"):
        self.data_directory = directory

    def set_user_directory(self, directory):
        self.user_directory = directory

    def set_data_file(self, filename):
        self.datafile = filename

    def startup(self):
        user_data_directory = self.use_folder
        self.set_data_directory(user_data_directory)
        print "Raw Data in :   ", user_data_directory
        user_data_processed = self.out_folder
        self.set_user_directory(directory=user_data_processed)
        print "Processed Data in :   ", user_data_processed

    # Returns the calibration filename
    def get_cal(self):
        return self.cal_dir + "WISH_cycle_10_3_noends_10to10.cal"
        # return "/home/mp43/Desktop/Si_Mar15/test_detOffsets_SiMar15_noends.cal"
        # return "/home/ryb18365/Desktop/WISH_cycle_10_3_noends_10to10_dodgytubesremoved.cal"

    # Returns the grouping filename
    def get_group_file(self):
        return self.cal_dir + "WISH_cycle_10_3_noends_10to10.cal"
        # return "/home/mp43/Desktop/Si_Mar15/test_detOffsets_SiMar15_noends.cal"
        # return "/home/ryb18365/Desktop/WISH_cycle_10_3_noends_10to10_dodgytubesremoved.cal"

    def get_vanadium(self, panel, cycle="09_4"):
        vanadium_string = {
                "09_2": "vana318-" + str(panel) + "foc-rmbins-smooth50.nx5",
                "09_3": "vana935-" + str(panel) + "foc-SS.nx5",
                "09_4": "vana3123-" + str(panel) + "foc-SS.nx5",
                "09_5": "vana3123-" + str(panel) + "foc-SS.nx5",
                "11_1": "vana17718-" + str(panel) + "foc-SS.nxs",
                "11_2": "vana16812-" + str(panel) + "foc-SS.nx5",
                "11_3": "vana18590-" + str(panel) + "foc-SS-new.nxs",
                "11_4": "vana38428-" + str(panel) + "foc-SF-SS.nxs",
                "18_2": "WISHvana41865-" + str(panel) + "foc.nxs"
        }
        return self.cal_dir + vanadium_string.get(cycle)

    def get_empty(self, panel, se="WISHcryo", cycle="09_4"):
        if se == "WISHcryo":
            empty_string = {
                "09_2": "emptycryo322-" + str(panel) + "-smooth50.nx5",
                "09_3": "emptycryo1725-" + str(panel) + "foc.nx5",
                "09_4": "emptycryo3307-" + str(panel) + "foc.nx5",
                "09_5": "emptycryo16759-" + str(panel) + "foc.nx5",
                "11_1": "emptycryo17712-" + str(panel) + "foc-SS.nxs",
                "11_2": "emptycryo16759-" + str(panel) + "foc-SS.nx5",
                "11_3": "emptycryo17712-" + str(panel) + "foc-SS-new.nxs",
                "11_4": "empty_mag20620-" + str(panel) + "foc-HR-SF.nxs"

            }
            return self.cal_dir + empty_string.get(cycle)

        if se == "candlestick":
            empty_string = {
                "09_3": "emptyinst1726-" + str(panel) + "foc-monitor.nxs",
                "09_4": "emptyinst3120-" + str(panel) + "foc.nxs",
                "11_4": "emptyinst19618-" + str(panel) + "foc-SF-S.nxs",
                "17_1": "emptyinst38581-" + str(panel) + "foc.nxs"
            }
            return self.cal_dir + empty_string.get(cycle)

    def generate_name_from_run(self, run_number, ext):
        filename = "WISH"+str(run_number)
        filename = filename+"."+ext
        print filename
        return filename

    def get_file_name(self, run_number, ext):
        if ext[0] != 's':
            data_dir = self.data_directory
        else:
            # data_dir="/datad/ndxwish/"
            data_dir = self.data_directory
        digit = len(str(run_number))
        filename = os.path.join(data_dir, "WISH")
        for i in range(0, 8 - digit):
            filename = filename + "0"
        filename += str(run_number) + "." + ext
        return filename

    def return_panel(self, panel):
        print panel
        return self.return_panel.get(panel)

    # Reads a wish data file return a workspace with a short name
    def read_file(self, number, panel, ext):
        if type(number) is int:
            filename = self.datafile  # Changed as full path is set in main now
            ext = filename.split('.')[-1]  # Get the extension from the inputted filename
            print "Extension is: " + ext
            # if (ext[0:10]=="nxs_event"):
            #    filename=WISH_getfilename(number,"nxs")
            print "will be reading filename..." + filename
            panel_min, panel_max = self.return_panel.get(panel)
            if panel != 0:
                output = "w" + str(number) + "_" + str(panel)
            else:
                output = "w" + str(number)
            shared_load_files(ext, filename, output, str(panel_max), str(panel_min), False)
            if ext == "nxs_event":
                mantid.LoadEventNexus(Filename=filename, OutputWorkspace=output, LoadMonitors='1')
                mantid.RenameWorkspace(output + "_monitors", "w" + str(number) + "_monitors")
                mantid.Rebin(InputWorkspace=output, OutputWorkspace=output, Params='6000,-0.00063,110000')
                panel_min, panel_max = self.return_panel.get(panel)
                mantid.CropWorkspace(InputWorkspace=output, OutputWorkspace=output, StartWorkspaceIndex=panel_min - 6,
                                     EndWorkspaceIndex=panel_max - 6)
                mantid.MaskBins(InputWorkspace=output,  OutputWorkspace=output, XMin=99900, XMax=106000)

                print "full nexus event file loaded"
            if ext[0:10] == "nxs_event_":
                label, tmin, tmax = split_string_event(ext)
                output = output + "_" + label
                if tmax == "end":
                    mantid.LoadEventNexus(Filename=filename, OutputWorkspace=output, FilterByTimeStart=tmin,
                                          LoadMonitors='1',
                                          MonitorsAsEvents='1', FilterMonByTimeStart=tmin)

                else:
                    mantid.LoadEventNexus(Filename=filename, OutputWorkspace=output, FilterByTimeStart=tmin,
                                          FilterByTimeStop=tmax,
                                          LoadMonitors='1', MonitorsAsEvents='1', FilterMonByTimeStart=tmin,
                                          FilterMonByTimeStop=tmax)
                    mantid.RenameWorkspace(output + "_monitors", "w" + str(number) + "_monitors")

                print "renaming monitors done!"
                mantid.Rebin(InputWorkspace=output, OutputWorkspace=output, Params='6000,-0.00063,110000')
                panel_min, panel_max = self.return_panel.get(panel)
                mantid.CropWorkspace(InputWorkspace=output, OutputWorkspace=output, StartWorkspaceIndex=panel_min - 6,
                                     EndWorkspaceIndex=panel_max - 6)
                mantid.MaskBins(InputWorkspace=output,  OutputWorkspace=output, XMin=99900, XMax=106000)
                print "nexus event file chopped"

        else:
            n1, n2 = split_string(number)
            output = "w" + str(n1) + "_" + str(n2) + "-" + str(panel)
            filename = self.get_file_name(n1, ext)
            print "reading filename..." + filename
            panel_min, panel_max = self.return_panel.get(panel)
            output1 = "w" + str(n1) + "-" + str(panel)
            mantid.LoadRaw(Filename=filename, OutputWorkspace=output1, SpectrumMin=str(panel_min),
                           SpectrumMax=str(panel_max), LoadLogFiles="0")
            filename = self.get_file_name(n2, ext)
            print "reading filename..." + filename
            panel_min, panel_max = self.return_panel.get(panel)
            output2 = "w" + str(n2) + "-" + str(panel)
            mantid.LoadRaw(Filename=filename, OutputWorkspace=output2, SpectrumMin=str(panel_min),
                           SpectrumMax=str(panel_max), LoadLogFiles="0")
            mantid.MergeRuns(InputWorkspaces=[output1, output2], OutputWorkspace=output)
            mantid.DeleteWorkspace(output1)
            mantid.DeleteWorkspace(output2)
            mantid.ConvertUnits(InputWorkspace=output, OutputWorkspace=output, Target="Wavelength", Emode="Elastic")

        lambda_min, lambda_max = get_lambda_range()
        mantid.CropWorkspace(InputWorkspace=output, OutputWorkspace=output, XMin=lambda_min, XMax=lambda_max)
        monitor = self.process_incident_monitor(number, ext)
        print "first norm to be done"
        mantid.NormaliseToMonitor(InputWorkspace=output, OutputWorkspace=output + "norm1", MonitorWorkspace=monitor)
        print "second norm to be done"
        mantid.NormaliseToMonitor(InputWorkspace=output + "norm1", OutputWorkspace=output + "norm2",
                                  MonitorWorkspace=monitor, IntegrationRangeMin=0.7, IntegrationRangeMax=10.35)
        mantid.DeleteWorkspace(output)
        mantid.DeleteWorkspace(output + "norm1")
        mantid.RenameWorkspace(InputWorkspace=output + "norm2", OutputWorkspace=output)
        mantid.ConvertUnits(InputWorkspace=output, OutputWorkspace=output, Target="TOF", EMode="Elastic")
        mantid.ReplaceSpecialValues(InputWorkspace=output, OutputWorkspace=output, NaNValue=0.0, NaNError=0.0,
                                    InfinityValue=0.0, InfinityError=0.0)
        return output

    # Focus data set for a given panel and return the workspace
    def focus_one_panel(self, work, focus, panel):
        panel_crop = {
            1: (0.8, 53.3),
            2: (0.5, 13.1),
            3: (0.5, 7.77),
            4: (0.4, 5.86),
            5: (0.35, 4.99),
            6: (0.35, 4.99),
            7: (0.4, 5.86),
            8: (0.5, 7.77),
            9: (0.5, 13.1),
            10: (0.8, 53.3)
        }
        d_min, d_max = panel_crop.get(panel)
        mantid.AlignAndFocusPowder(InputWorkspace=work, OutputWorkspace=focus, GroupFileName=self.get_group_file(),
                                   CalFileName=self.get_cal(), RunRebin=False)
        mantid.ConvertUnits(InputWorkspace=focus, OutputWorkspace=focus, Target="dSpacing", EMode="Elastic")
        mantid.CropWorkspace(InputWorkspace=focus, OutputWorkspace=focus, XMin=d_min, XMax=d_max)
        return focus

    def focus(self, work, panel):
        focus = work + "foc"
        if panel != 0:
            self.focus_one_panel(work, focus, panel)
        else:
            self.focus_one_panel(work, focus, panel)
            split_workspace(focus)

    def process(self, number, panel, ext, se_sample="WISHcryo", empty_se_cycle="09_4", cycle_vanadium="09_4",
                absorb=False, nd=0.0, xs=0.0, xa=0.0, h=0.0, r=0.0):
        w = self.read_file(number, panel, ext)
        print "file read and normalized"
        self.focus(w, panel)
        print "focussing done!"
        if absorb:
            print absorb
            mantid.ConvertUnits(InputWorkspace=w, OutputWorkspace=w, Target="Wavelength", EMode="Elastic")
            mantid.CylinderAbsorption(InputWorkspace=w, OutputWorkspace="T",
                                      CylinderSampleHeight=h, CylinderSampleRadius=r, AttenuationXSection=xa,
                                      ScatteringXSection=xs, SampleNumberDensity=nd,
                                      NumberOfSlices="10", NumberOfAnnuli="10", NumberOfWavelengthPoints="25",
                                      ExpMethod="Normal")
            mantid.Divide(LHSWorkspace=w, RHSWorkspace="T", OutputWorkspace=w)
            mantid.DeleteWorkspace("T")
            mantid.ConvertUnits(InputWorkspace=w, OutputWorkspace=w, Target="dSpacing", EMode="Elastic")
            print "absorb done"

        if type(number) is int:
            focused_workspace_name = "w" + str(number) + "_" + str(panel) + "foc"
            if len(ext) > 9:
                label, tmin, tmax = split_string_event(ext)
                focused_workspace_name = "w" + str(number) + "_" + str(panel) + "_" + label + "foc"
        else:
            n1, n2 = split_string(number)
            focused_workspace_name = "w" + str(n1) + "_" + str(n2) + "_" + str(panel) + "foc"
        if self.deleteWorkspace:
            mantid.DeleteWorkspace(w)
        print panel
        print se_sample
        print empty_se_cycle
        print self.get_empty(panel, se_sample, empty_se_cycle)
        if panel == 0:
            for i in range(1, 11):
                focused_workspace_name = "w" + str(number) + "_" + str(i) + "foc"
                mantid.LoadNexusProcessed(Filename=self.get_empty(i, se_sample, empty_se_cycle),
                                          OutputWorkspace="empty")
                mantid.RebinToWorkspace(WorkspaceToRebin="empty", WorkspaceToMatch=focused_workspace_name,
                                        OutputWorkspace="empty")
                mantid.Minus(LHSWorkspace=focused_workspace_name, RHSWorkspace="empty",
                             OutputWorkspace=focused_workspace_name)
                mantid.DeleteWorkspace("empty")
                print "will try to load a vanadium with the name:" + self.get_vanadium(i, cycle_vanadium)
                mantid.LoadNexusProcessed(Filename=self.get_vanadium(i, cycle_vanadium), OutputWorkspace="vana")
                mantid.RebinToWorkspace(WorkspaceToRebin="vana", WorkspaceToMatch=focused_workspace_name,
                                        OutputWorkspace="vana")
                mantid.Divide(LHSWorkspace=focused_workspace_name, RHSWorkspace="vana",
                              OutputWorkspace=focused_workspace_name)
                mantid.DeleteWorkspace("vana")
                mantid.ConvertUnits(InputWorkspace=focused_workspace_name, OutputWorkspace=focused_workspace_name,
                                    Target="TOF", EMode="Elastic")
                mantid.ReplaceSpecialValues(InputWorkspace=focused_workspace_name,
                                            OutputWorkspace=focused_workspace_name, NaNValue=0.0,
                                            NaNError=0.0,
                                            InfinityValue=0.0, InfinityError=0.0)
                mantid.SaveGSS(InputWorkspace=focused_workspace_name,
                               Filename=os.path.join(self.user_directory, (str(number) + "-" + str(i) + ext
                                                                                + ".gss")),
                               Append=False, Bank=1)
                mantid.SaveFocusedXYE(focused_workspace_name,
                                      os.path.join(self.user_directory, (str(number) + "-" + str(i) + ext
                                                                              + ".dat")))
                mantid.SaveNexusProcessed(focused_workspace_name,
                                          os.path.join(self.user_directory, (str(number) + "-" + str(i) + ext
                                                                                  + ".nxs")))
        else:
            mantid.LoadNexusProcessed(Filename=self.get_empty(panel, se_sample, empty_se_cycle),
                                      OutputWorkspace="empty")
            mantid.RebinToWorkspace(WorkspaceToRebin="empty", WorkspaceToMatch=focused_workspace_name,
                                    OutputWorkspace="empty")
            mantid.Minus(LHSWorkspace=focused_workspace_name, RHSWorkspace="empty",
                         OutputWorkspace=focused_workspace_name)
            mantid.DeleteWorkspace("empty")
            print "will try to load a vanadium with the name:" + self.get_vanadium(panel, cycle_vanadium)
            mantid.LoadNexusProcessed(Filename=self.get_vanadium(panel, cycle_vanadium), OutputWorkspace="vana")
            mantid.RebinToWorkspace(WorkspaceToRebin="vana", WorkspaceToMatch=focused_workspace_name,
                                    OutputWorkspace="vana")
            mantid.Divide(LHSWorkspace=focused_workspace_name, RHSWorkspace="vana",
                          OutputWorkspace=focused_workspace_name)
            mantid.DeleteWorkspace("vana")
            mantid.ConvertUnits(InputWorkspace=focused_workspace_name, OutputWorkspace=focused_workspace_name,
                                Target="TOF", EMode="Elastic")
            mantid.ReplaceSpecialValues(InputWorkspace=focused_workspace_name, OutputWorkspace=focused_workspace_name,
                                        NaNValue=0.0, NaNError=0.0,
                                        InfinityValue=0.0, InfinityError=0.0)
            mantid.SaveGSS(InputWorkspace=focused_workspace_name,
                           Filename=os.path.join(self.user_directory, (str(number) + "-" + str(panel) + ext
                                                 + ".gss")),
                           Append=False, Bank=1)
            mantid.SaveFocusedXYE(focused_workspace_name,
                                  os.path.join(self.user_directory, (str(number) + "-" + str(panel) + ext
                                               + ".dat")))
            mantid.SaveNexusProcessed(focused_workspace_name,
                                      os.path.join(self.user_directory, (str(number) + "-" + str(panel) + ext
                                                   + ".nxs")))
        return focused_workspace_name

    # Create a corrected vanadium (normalise,corrected for attenuation and empty, strip peaks) and
    # save a a nexus processed file.
    # It looks like smoothing of 100 works quite well
    def create_normalised_vanadium(self, van, empty, panel, vh, vr, cycle_van="18_2", cycle_empty="17_1"):
        self.startup()
        self.set_data_file(self.generate_name_from_run(van, "nxs"))
        self.set_data_directory("/archive/ndxwish/Instrument/data/cycle_" + cycle_van + "/")
        van = self.read_file(van, panel, "nxs_event")
        self.startup()
        self.set_data_file(self.generate_name_from_run(empty, "nxs"))
        self.set_data_directory("/archive/ndxwish/Instrument/data/cycle_" + cycle_empty + "/")
        empty = self.read_file(empty, panel, "nxs_event")
        mantid.Minus(LHSWorkspace=van, RHSWorkspace=empty, OutputWorkspace=van)
        print "read van and empty"
        mantid.DeleteWorkspace(empty)
        mantid.ConvertUnits(InputWorkspace=van, OutputWorkspace=van, Target="Wavelength", EMode="Elastic")
        mantid.CylinderAbsorption(InputWorkspace=van, OutputWorkspace="T",
                                  CylinderSampleHeight=str(vh), CylinderSampleRadius=str(vr),
                                  AttenuationXSection="4.8756",
                                  ScatteringXSection="5.16", SampleNumberDensity="0.07118",
                                  NumberOfSlices="10", NumberOfAnnuli="10", NumberOfWavelengthPoints="25",
                                  ExpMethod="Normal")
        mantid.Divide(LHSWorkspace=van, RHSWorkspace="T", OutputWorkspace=van)
        mantid.DeleteWorkspace("T")
        mantid.ConvertUnits(InputWorkspace=van, OutputWorkspace=van, Target="TOF", EMode="Elastic")
        self.focus(van, panel)
        mantid.DeleteWorkspace(van)
        return

    def create_empty(self, empty, panel):
        empty = self.read_file(empty, panel, "raw")
        self.focus(empty, panel)
        return

    # Have made no changes here as not called (may not work in future though)
    def monitors(self, rb, ext):
        # data_dir = WISH_dir()
        filename = self.get_file_name(rb, ext)
        workspace_out = "w" + str(rb)
        print "reading File..." + filename
        mantid.LoadRaw(Filename=filename, OutputWorkspace=workspace_out, SpectrumMin=str(1), SpectrumMax=str(5),
                       LoadLogFiles="0")
        mantid.NormaliseByCurrent(InputWorkspace=workspace_out, OutputWorkspace=workspace_out)
        mantid.ConvertToDistribution(workspace_out)
        return workspace_out

    def process_incident_monitor(self, number, ext):
        print("process monitor")
        if type(number) is int:
            filename = self.datafile
            works = "monitor" + str(number)
            if shared_load_files(ext, filename, works, 4, 4, True):
                works = "monitor" + str(number)
            if ext[0:9] == "nxs_event":
                temp = "w" + str(number) + "_monitors"
                works = "w" + str(number) + "_monitor4"
                mantid.Rebin(InputWorkspace=temp, OutputWorkspace=temp, Params='6000,-0.00063,110000',
                             PreserveEvents=False)
                mantid.ExtractSingleSpectrum(InputWorkspace=temp, OutputWorkspace=works, WorkspaceIndex=3)
        else:
            n1, n2 = split_string(number)
            works = "monitor" + str(n1) + "_" + str(n2)
            filename = self.get_file_name(n1, ext)
            works1 = "monitor" + str(n1)
            mantid.LoadRaw(Filename=filename, OutputWorkspace=works1, SpectrumMin=4, SpectrumMax=4, LoadLogFiles="0")
            filename = self.get_file_name(n2, ext)
            works2 = "monitor" + str(n2)
            mantid.LoadRaw(Filename=filename, OutputWorkspace=works2, SpectrumMin=4, SpectrumMax=4, LoadLogFiles="0")
            mantid.MergeRuns(InputWorkspaces=[works1, works2], OutputWorkspace=works)
            mantid.DeleteWorkspace(works1)
            mantid.DeleteWorkspace(works2)
            mantid.ConvertUnits(InputWorkspace=works, OutputWorkspace=works, Target="Wavelength", Emode="Elastic")
        l_min, l_max = get_lambda_range()
        mantid.CropWorkspace(InputWorkspace=works, OutputWorkspace=works, XMin=l_min, XMax=l_max)
        ex_regions = n.array([[4.57, 4.76],
                              [3.87, 4.12],
                              [2.75, 2.91],
                              [2.24, 2.50]])
        mantid.ConvertToDistribution(works)
        for reg in range(0, 4):
            mantid.MaskBins(InputWorkspace=works, OutputWorkspace=works, XMin=ex_regions[reg, 0],
                            XMax=ex_regions[reg, 1])
        mantid.ConvertFromDistribution(works)
        return works

    def minus_empty_cans(self, run_number, empty):
        panel_list = ['-1foc', '-2foc', '-3foc', '-4foc', '-5foc', '-6foc', '-7foc', '-8foc', '-9foc', '-10foc',
                      '-1_10foc',
                      '-2_9foc', '-3_8foc', '-4_7foc', '-5_6foc']
        for p in panel_list:
            mantid.Minus(LHSWorkspace='w' + str(run_number) + p, RHSWorkspace='w' + str(empty) + p,
                         OutputWorkspace='w' + str(run_number) + 'minus' + str(empty) + p)
            mantid.ConvertUnits(InputWorkspace='w' + str(run_number) + 'minus' + str(empty) + p,
                                OutputWorkspace='w' + str(run_number) + 'minus' + str(empty) + p
                                                + '-d', Target='dSpacing')
            mantid.SaveGSS("w" + str(run_number) + 'minus' + str(empty) + p,
                           os.path.join(self.user_directory, (str(run_number) + p + ".gss")), Append=False, Bank=1)
        return

    def main(self):
        self.validate()
        print(self.user_directory)
        print(self.datafile)
        i = get_run_number(self.datafile)
        for panel in range(1, 11):
            output_workspace = self.process(i, panel, "raw", "candlestick", "17_1", "18_2", absorb=False, nd=0.0,
                                            xs=0.0, xa=0.0, h=4.0, r=0.4)
            mantid.ConvertUnits(InputWorkspace=output_workspace, OutputWorkspace=output_workspace + "-d",
                                Target="dSpacing", EMode="Elastic")
        for panel in range(1, 6):
            self.save_combined_panel(i, panel)

    def save_combined_panel(self, run, panel):
        panel_combination = {
            5: "6",
            4: "7",
            3: "8",
            2: "9",
            1: "10"
        }
        combined_panel = str(panel) + "-" + panel_combination.get(panel)
        input_workspace1 = "w" + str(run) + "_" + str(panel) + "foc"
        input_workspace2 = "w" + str(run) + "_" + panel_combination.get(panel) + "foc"
        combined_workspace = "w" + str(run) + "_" + combined_panel + "foc"
        mantid.RebinToWorkspace(WorkspaceToRebin=input_workspace2, WorkspaceToMatch=input_workspace1,
                                OutputWorkspace=input_workspace2, PreserveEvents='0')
        mantid.Plus(LHSWorkspace=input_workspace1, RHSWorkspace=input_workspace2, OutputWorkspace=combined_workspace)
        mantid.ConvertUnits(InputWorkspace=combined_workspace, OutputWorkspace=combined_workspace + "-d",
                            Target="dSpacing", EMode="Elastic")
        mantid.SaveGSS(combined_workspace, os.path.join(self.user_directory,
                                                        (str(run) + "-" + combined_panel + "raw" + ".gss")),
                       Append=False, Bank=1)
        mantid.SaveFocusedXYE(combined_workspace, os.path.join(self.user_directory,
                                                               (str(run) + "-" + combined_panel + "raw" + ".dat")))
        mantid.SaveNexusProcessed(combined_workspace, os.path.join(self.user_directory,
                                                                   (str(run) + "-" + combined_panel + "raw" + ".nxs")))

    # test vanadium is 41870 test empty is 38581
    def create_vanadium(self, vanadium_run, empty_run):
        # ######### use the lines below to process a LoadRawvanadium run                               #################
        for j in range(1, 11):
            self.create_normalised_vanadium(vanadium_run, empty_run, j, 4.0, 0.15, cycle_van="18_2", cycle_empty="17_1")
            vanadium_workspace = str(vanadium_run)+"_"+str(j)+"foc"
            mantid.CropWorkspace(InputWorkspace="w"+vanadium_workspace, OutputWorkspace="w"+vanadium_workspace,
                                 XMin='0.35', XMax='5.0')
            remove_peaks_spline_smooth_empty("w"+vanadium_workspace, j)
            mantid.SaveNexusProcessed("w"+vanadium_workspace,
                                      os.path.join(self.user_directory, ("vana41865-" + str(j) + "foc.nxs")))#"vana" + vanadium_workspace + ".nxs")))

    def run_script(self, run):
        if self.name == "__main__":
            self.startup()
            self.set_data_file(self.get_file_name(run, "nxs"))

            self.main()


def split_string(t):
    index_p = 0
    for i in range(0, len(t)):
        if t[i] == "+":
            index_p = i
    if index_p != 0:
        return int(t[0:index_p]), int(t[index_p + 1:len(t)])


def split_string_event(t):
    # this assumes the form nxs_event_label_tmin_tmax
    index_ = []
    for i in range(10, len(t)):
        if t[i] == "_":
            index_.append(i)
    label = t[10:index_[0]]
    t_min = t[index_[0] + 1:index_[1]]
    t_max = t[index_[1] + 1:len(t)]
    return label, t_min, t_max


def split_workspace(focus):
    for i in range(0, 11):
        out = focus[0:len(focus) - 3] + "-" + str(i + 1) + "foc"
        mantid.ExtractSingleSpectrum(InputWorkspace=focus, OutputWorkspace=out, WorkspaceIndex=i)
        mantid.DeleteWorkspace(focus)
    return


# Get the run number from the input data path
def get_run_number(data_path):
    data_path = data_path.split('/')[-1]  # Get the file name
    data_path = data_path.split('.')[0]  # Remove the extension
    data_path = data_path[4:]  # Remove the WISH prefix
    return int(data_path)


# Get the valid wavelength range, i.e. excluding regions where choppers cut
def get_lambda_range():
    return 0.7, 10.35


def remove_peaks_spline_smooth_empty(works, panel):
    smooth_terms = {
        1: 30,
        2: 10,
        3: 15,
        4: 15,
        5: 10,
        10: 30,
        9: 10,
        8: 15,
        7: 15,
        6: 10
    }
    mantid.SmoothData(InputWorkspace=works, OutputWorkspace=works, NPoints=smooth_terms.get(panel))
    return works


def shared_load_files(ext, filename, output, spectrum_max, spectrum_min, is_monitor):
    if ext == "nxs":
        mantid.LoadNexus(Filename=filename, OutputWorkspace=output, SpectrumMin=spectrum_min,
                         SpectrumMax=spectrum_max)
        mantid.Rebin(InputWorkspace=output, OutputWorkspace=output, Params='6000,-0.00063,110000')
        if not is_monitor:
            mantid.MaskBins(InputWorkspace=output,  OutputWorkspace=output, XMin=99900, XMax=106000)
        mantid.ConvertUnits(InputWorkspace=output, OutputWorkspace=output, Target="Wavelength", Emode="Elastic")

        return True
    if ext == "raw":
        mantid.LoadRaw(Filename=filename, OutputWorkspace=output, SpectrumMin=spectrum_min,
                       SpectrumMax=spectrum_max, LoadLogFiles="0")
        mantid.Rebin(InputWorkspace=output, OutputWorkspace=output, Params='6000,-0.00063,110000')
        if not is_monitor:
            mantid.MaskBins(InputWorkspace=output,  OutputWorkspace=output, XMin=99900, XMax=106000)
        return True
    if ext[0] == "s":
        mantid.LoadRaw(Filename=filename, OutputWorkspace=output, SpectrumMin=spectrum_min,
                       SpectrumMax=spectrum_max, LoadLogFiles="0")
        mantid.Rebin(InputWorkspace=output, OutputWorkspace=output, Params='6000,-0.00063,110000')
        if not is_monitor:
            mantid.MaskBins(InputWorkspace=output,  OutputWorkspace=output, XMin=99900, XMax=106000)
        return True
    return False
