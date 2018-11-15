import numpy as n
import os
import mantid.simpleapi as mantid


class Wish:
    NUM_PANELS = 11
    NUM_MONITORS = 6

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

    # Returns the grouping filename
    def get_group_file(self):
        return self.cal_dir + "WISH_cycle_10_3_noends_10to10.cal"

    def get_vanadium(self, panel, cycle="09_4"):
        vanadium_string = {
            "09_2": "vana318-{}foc-rmbins-smooth50.nx5",
            "09_3": "vana935-{}foc-SS.nx5",
            "09_4": "vana3123-{}foc-SS.nx5",
            "09_5": "vana3123-{}foc-SS.nx5",
            "11_1": "vana17718-{}foc-SS.nxs",
            "11_2": "vana16812-{}foc-SS.nx5",
            "11_3": "vana18590-{}foc-SS-new.nxs",
            "11_4": "vana38428-{}foc-SF-SS.nxs",
            "18_2": "WISHvana41865-{}foc.nxs"
        }
        return self.cal_dir + vanadium_string.get(cycle).format(panel)

    def get_empty(self, panel, se="WISHcryo", cycle="09_4"):
        if se == "WISHcryo":
            empty_string = {
                "09_2": "emptycryo322-{}-smooth50.nx5",
                "09_3": "emptycryo1725-{}foc.nx5",
                "09_4": "emptycryo3307-{}foc.nx5",
                "09_5": "emptycryo16759-{}foc.nx5",
                "11_1": "emptycryo17712-{}foc-SS.nxs",
                "11_2": "emptycryo16759-{}foc-SS.nx5",
                "11_3": "emptycryo17712-{}foc-SS-new.nxs",
                "11_4": "empty_mag20620-{}foc-HR-SF.nxs"

            }
            return self.cal_dir + empty_string.get(cycle).format(panel)

        if se == "candlestick":
            empty_string = {
                "09_3": "emptyinst1726-{}foc-monitor.nxs",
                "09_4": "emptyinst3120-{}foc.nxs",
                "11_4": "emptyinst19618-{}foc-SF-S.nxs",
                "17_1": "emptyinst38581-{}foc.nxs"
            }
            return self.cal_dir + empty_string.get(cycle).format(panel)

    def get_file_name(self, run_number, extension):
        if extension[0] != 's':
            data_dir = self.data_directory
        else:
            data_dir = self.data_directory
        digit = len(str(run_number))
        filename = os.path.join(data_dir, "WISH")
        for i in range(8 - digit):
            filename = filename + "0"
        filename += str(run_number) + "." + extension
        return filename

    def return_panel(self, panel):
        print panel
        return self.return_panel.get(panel)

    # Reads a wish data file return a workspace with a short name
    def read_file(self, number, panel, extension):
        if type(number) is int:
            filename = self.datafile  # Changed as full path is set in main now
            extension = filename.split('.')[-1]  # Get the extension from the inputted filename
            print "Extension is: " + extension
            print "will be reading filename..." + filename
            panel_min, panel_max = self.return_panel.get(panel)
            if panel != 0:
                output = "w{0}_{1}".format(number, panel)
            else:
                output = "w{}".format(number)
            shared_load_files(extension, filename, output, str(panel_max), str(panel_min), False)
            if extension == "nxs_event":
                mantid.LoadEventNexus(Filename=filename, OutputWorkspace=output, LoadMonitors='1')
                mantid.RenameWorkspace(output + "_monitors", "w{}_monitors".format(number))
                mantid.Rebin(InputWorkspace=output, OutputWorkspace=output, Params='6000,-0.00063,110000')
                panel_min, panel_max = self.return_panel.get(panel)
                mantid.CropWorkspace(InputWorkspace=output, OutputWorkspace=output, StartWorkspaceIndex=panel_min - 6,
                                     EndWorkspaceIndex=panel_max - 6)
                mantid.MaskBins(InputWorkspace=output, OutputWorkspace=output, XMin=99900, XMax=106000)

                print "full nexus event file loaded"
            if extension[:10] == "nxs_event_":
                label, time_min, time_max = split_string_event(extension)
                output = output + "_" + label
                if time_max == "end":
                    mantid.LoadEventNexus(Filename=filename, OutputWorkspace=output, FilterByTimeStart=time_min,
                                          LoadMonitors='1',
                                          MonitorsAsEvents='1', FilterMonByTimeStart=time_min)

                else:
                    mantid.LoadEventNexus(Filename=filename, OutputWorkspace=output, FilterByTimeStart=time_min,
                                          FilterByTimeStop=time_max,
                                          LoadMonitors='1', MonitorsAsEvents='1', FilterMonByTimeStart=time_min,
                                          FilterMonByTimeStop=time_max)
                    mantid.RenameWorkspace(output + "_monitors", "w{}_monitors".format(number))

                print "renaming monitors done!"
                mantid.Rebin(InputWorkspace=output, OutputWorkspace=output, Params='6000,-0.00063,110000')
                panel_min, panel_max = self.return_panel.get(panel)
                mantid.CropWorkspace(InputWorkspace=output, OutputWorkspace=output, StartWorkspaceIndex=panel_min
                                     - Wish.NUM_MONITORS, EndWorkspaceIndex=panel_max - Wish.NUM_MONITORS)
                mantid.MaskBins(InputWorkspace=output, OutputWorkspace=output, XMin=99900, XMax=106000)

                print "nexus event file chopped"

        else:
            first_number, second_number = split_string(number)
            output = "w{0}_{1}-{2}".format(first_number, second_number, panel)
            filename = self.get_file_name(first_number, extension)
            print "reading filename..." + filename
            panel_min, panel_max = self.return_panel.get(panel)
            output1 = "w{0}-{1}".format(first_number, panel)
            mantid.LoadRaw(Filename=filename, OutputWorkspace=output1, SpectrumMin=str(panel_min),
                           SpectrumMax=str(panel_max), LoadLogFiles="0")
            filename = self.get_file_name(second_number, extension)
            print "reading filename..." + filename
            panel_min, panel_max = self.return_panel.get(panel)
            output2 = "w{0}-{1}".format(second_number, panel)
            mantid.LoadRaw(Filename=filename, OutputWorkspace=output2, SpectrumMin=str(panel_min),
                           SpectrumMax=str(panel_max), LoadLogFiles="0")
            mantid.MergeRuns(InputWorkspaces=[output1, output2], OutputWorkspace=output)
            mantid.DeleteWorkspace(output1)
            mantid.DeleteWorkspace(output2)
            mantid.ConvertUnits(InputWorkspace=output, OutputWorkspace=output, Target="Wavelength", Emode="Elastic")

        lambda_min, lambda_max = get_lambda_range()
        mantid.CropWorkspace(InputWorkspace=output, OutputWorkspace=output, XMin=lambda_min, XMax=lambda_max)
        monitor = self.process_incident_monitor(number, extension)
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
        mantid.CropWorkspaceRagged(InputWorkspace=focus, OutputWorkspace=focus, XMin=d_min, XMax=d_max)
        return focus

    def focus(self, work, panel):
        focus = work + "foc"
        if panel != 0:
            self.focus_one_panel(work, focus, panel)
        else:
            self.focus_one_panel(work, focus, panel)
            split_workspace(focus)

    def process(self, number, panel, extension, se_sample="WISHcryo", empty_se_cycle="09_4", cycle_vanadium="09_4",
                absorb=False, number_density=0.0, scattering_x=0.0, attenuation_x=0.0, cylinder_height=0.0,
                cylinder_radius=0.0):
        workspace_to_focus = self.read_file(number, panel, extension)
        print "file read and normalized"
        if absorb:
            absorption_corrections(height=cylinder_height, number_density=number_density, radius=cylinder_radius,
                                   input_workspace=workspace_to_focus, attenuation_x=attenuation_x,
                                   scattering_x=scattering_x)
            mantid.ConvertUnits(InputWorkspace=workspace_to_focus, OutputWorkspace=workspace_to_focus, Target="TOF",
                                EMode="Elastic")
        self.focus(workspace_to_focus, panel)
        print "focussing done!"

        if type(number) is int:
            focused_workspace_name = "w{0}_{1}foc".format(number, panel)
            if len(extension) > 9:
                label, tmin, tmax = split_string_event(extension)
                focused_workspace_name = "w{0}_{1}_{2}foc".format(number, panel, label)
        else:
            first_number, second_number = split_string(number)
            focused_workspace_name = "w{0}_{1}_{2}foc".format(first_number, second_number, panel)
        if self.deleteWorkspace:
            mantid.DeleteWorkspace(workspace_to_focus)
        if panel == 0:
            for i in range(1, Wish.NUM_PANELS):
                focused_workspace_name = "w{0}_{1}foc".format(number, i)
                self.empty_and_van_correction(cycle_vanadium, empty_se_cycle, extension, focused_workspace_name, number,
                                              i, se_sample)
        else:
            self.empty_and_van_correction(cycle_vanadium, empty_se_cycle, extension, focused_workspace_name, number,
                                          panel, se_sample)
        return focused_workspace_name

    def empty_and_van_correction(self, cycle_vanadium, empty_se_cycle, extension, focused_workspace_name, number, panel,
                                 se_sample):

        mantid.LoadNexusProcessed(Filename=self.get_empty(panel, se_sample, empty_se_cycle),
                                  OutputWorkspace="empty")
        mantid.RebinToWorkspace(WorkspaceToRebin="empty", WorkspaceToMatch=focused_workspace_name,
                                OutputWorkspace="empty")
        mantid.Minus(LHSWorkspace=focused_workspace_name, RHSWorkspace="empty", OutputWorkspace=focused_workspace_name)
        mantid.DeleteWorkspace("empty")

        print "will try to load a vanadium with the name:" + self.get_vanadium(panel, cycle_vanadium)
        mantid.LoadNexusProcessed(Filename=self.get_vanadium(panel, cycle_vanadium), OutputWorkspace="vana")
        mantid.RebinToWorkspace(WorkspaceToRebin="vana", WorkspaceToMatch=focused_workspace_name,
                                OutputWorkspace="vana")
        mantid.Divide(LHSWorkspace=focused_workspace_name, RHSWorkspace="vana", OutputWorkspace=focused_workspace_name)
        mantid.DeleteWorkspace("vana")
        mantid.ConvertUnits(InputWorkspace=focused_workspace_name, OutputWorkspace=focused_workspace_name, Target="TOF",
                            EMode="Elastic")
        mantid.ReplaceSpecialValues(InputWorkspace=focused_workspace_name, OutputWorkspace=focused_workspace_name,
                                    NaNValue=0.0, NaNError=0.0, InfinityValue=0.0, InfinityError=0.0)

        mantid.SaveGSS(InputWorkspace=focused_workspace_name,
                       Filename=os.path.join(self.user_directory, ("{0}-{1}{2}.gss".format(number, panel,  extension))),
                       Append=False, Bank=1)
        mantid.SaveFocusedXYE(focused_workspace_name,
                              os.path.join(self.user_directory, (str(number) + "-" + str(panel) + extension + ".dat")))
        mantid.SaveNexusProcessed(focused_workspace_name,
                                  os.path.join(self.user_directory,
                                               (str(number) + "-" + str(panel) + extension + ".nxs")))

    # Create a corrected vanadium (normalise,corrected for attenuation and empty, strip peaks) and
    # save a a nexus processed file.
    # It looks like smoothing of 100 works quite well
    def create_normalised_vanadium(self, van, empty, panel, vanadium_height, vanadium_radius, cycle_van="18_2",
                                   cycle_empty="17_1"):
        self.startup()
        self.set_data_file(generate_name_from_run(van, "nxs"))
        directories = "/archive/ndxwish/Instrument/data/cycle_{0}/"
        self.set_data_directory(directories.format(cycle_van))
        van = self.read_file(van, panel, "nxs_event")
        self.startup()
        self.set_data_file(generate_name_from_run(empty, "nxs"))
        self.set_data_directory(directories.format(cycle_empty))
        empty = self.read_file(empty, panel, "nxs_event")
        mantid.Minus(LHSWorkspace=van, RHSWorkspace=empty, OutputWorkspace=van)
        print "read van and empty"
        mantid.DeleteWorkspace(empty)
        absorption_corrections(height=vanadium_height, number_density=0.07118, radius=vanadium_radius,
                               input_workspace=van, attenuation_x=4.8756, scattering_x=5.16)
        mantid.ConvertUnits(InputWorkspace=van, OutputWorkspace=van, Target="TOF", EMode="Elastic")
        self.focus(van, panel)
        mantid.DeleteWorkspace(van)
        return

    def create_empty(self, empty, panel):
        empty = self.read_file(empty, panel, "raw")
        self.focus(empty, panel)
        return

    def process_incident_monitor(self, number, extension):
        print("process monitor")
        if type(number) is int:
            filename = self.datafile
            works = "monitor" + str(number)
            if shared_load_files(extension, filename, works, 4, 4, True):
                works = "monitor" + str(number)
            if extension[:9] == "nxs_event":
                temp = "w{}_monitors".format(number)
                works = "w{}_monitor4".format(number)
                mantid.Rebin(InputWorkspace=temp, OutputWorkspace=temp, Params='6000,-0.00063,110000',
                             PreserveEvents=False)
                mantid.ExtractSingleSpectrum(InputWorkspace=temp, OutputWorkspace=works, WorkspaceIndex=3)
        else:
            n1, n2 = split_string(number)
            works = "monitor{0}_{1}".format(n1, n2)
            filename = self.get_file_name(n1, extension)
            works1 = "monitor{}".format(n1)
            mantid.LoadRaw(Filename=filename, OutputWorkspace=works1, SpectrumMin=4, SpectrumMax=4, LoadLogFiles="0")
            filename = self.get_file_name(n2, extension)
            works2 = "monitor{}".format(n2)
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
        for reg in range(4):
            mantid.MaskBins(InputWorkspace=works, OutputWorkspace=works, XMin=ex_regions[reg, 0],
                            XMax=ex_regions[reg, 1])
        mantid.ConvertFromDistribution(works)
        return works

    def minus_empty_cans(self, run_number, empty):
        panel_list = ['-1foc', '-2foc', '-3foc', '-4foc', '-5foc', '-6foc', '-7foc', '-8foc', '-9foc', '-10foc',
                      '-1_10foc', '-2_9foc', '-3_8foc', '-4_7foc', '-5_6foc']
        for p in panel_list:
            output = 'w{0}minus{1}{2}'.format(run_number, empty, p)
            mantid.Minus(LHSWorkspace='w{0}{1}'.format(run_number, p), RHSWorkspace='w{0}{1}'.format(empty, p),
                         OutputWorkspace=output)
            mantid.ConvertUnits(InputWorkspace=output, OutputWorkspace=output + '-d', Target='dSpacing')
            mantid.SaveGSS(output, os.path.join(self.user_directory, (str(run_number) + p + ".gss")), Append=False,
                           Bank=1)
        return

    def main(self):
        self.validate()
        print(self.user_directory)
        print(self.datafile)
        i = get_run_number(self.datafile)
        for panel in range(1, Wish.NUM_PANELS):
            output_workspace = self.process(i, panel, "raw", "candlestick", "17_1", "18_2", absorb=False,
                                            number_density=0.0, scattering_x=0.0, attenuation_x=0.0,
                                            cylinder_height=4.0, cylinder_radius=0.4)
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
        input_workspace1 = "w{0}_{1}foc".format(run, panel)
        input_workspace2 = "w{0}_{1}foc".format(run, panel_combination.get(panel))
        combined = "{0}{1}_{2}-{3}{4}".format("{0}", run, panel, panel_combination.get(panel), "{1}")
        combined_save = combined.format("", "{}")
        combined_ws = combined.format("w", "")

        mantid.RebinToWorkspace(WorkspaceToRebin=input_workspace2, WorkspaceToMatch=input_workspace1,
                                OutputWorkspace=input_workspace2, PreserveEvents='0')
        mantid.Plus(LHSWorkspace=input_workspace1, RHSWorkspace=input_workspace2, OutputWorkspace=combined_ws)
        mantid.ConvertUnits(InputWorkspace=combined_ws, OutputWorkspace=combined_ws+"-d", Target="dSpacing",
                            EMode="Elastic")

        mantid.SaveGSS(combined_ws, os.path.join(self.user_directory, combined_save.format("raw.gss")), Append=False,
                       Bank=1)
        mantid.SaveFocusedXYE(combined_ws, os.path.join(self.user_directory, combined_save.format("raw.dat")))
        mantid.SaveNexusProcessed(combined_ws, os.path.join(self.user_directory, combined_save.format("raw.nxs")))

    # test vanadium is 41870 test empty is 38581
    def create_vanadium(self, vanadium_run, empty_run):
        # ######### use the lines below to process a LoadRawvanadium run                               #################
        for panel in range(1, Wish.NUM_PANELS):
            self.create_normalised_vanadium(vanadium_run, empty_run, panel, 4.0, 0.15, cycle_van="18_2",
                                            cycle_empty="17_1")
            vanadium = "{0}{1}_{2}foc{3}"
            vanadium_workspace = vanadium.format("w", vanadium_run, panel, "")
            vanadium_save = vanadium.format("vana", vanadium_run, panel, ".nxs")
            mantid.CropWorkspace(InputWorkspace=vanadium_workspace, OutputWorkspace=vanadium_workspace,
                                 XMin='0.35', XMax='5.0')
            remove_peaks_spline_smooth_empty("w" + vanadium_workspace, panel)
            mantid.SaveNexusProcessed(vanadium_workspace,
                                      os.path.join(self.user_directory, vanadium_save))

    def run_script(self, run):
        if self.name == "__main__":
            self.startup()
            self.set_data_file(self.get_file_name(run, "nxs"))

            self.main()


def absorption_corrections(height, number_density, radius, input_workspace, attenuation_x, scattering_x):
    mantid.ConvertUnits(InputWorkspace=input_workspace, OutputWorkspace=input_workspace, Target="Wavelength",
                        EMode="Elastic")
    mantid.CylinderAbsorption(InputWorkspace=input_workspace, OutputWorkspace="T",
                              CylinderSampleHeight=height, CylinderSampleRadius=radius,
                              AttenuationXSection=attenuation_x,
                              ScatteringXSection=scattering_x, SampleNumberDensity=number_density,
                              NumberOfSlices="10", NumberOfAnnuli="10", NumberOfWavelengthPoints="25",
                              ExpMethod="Normal")
    mantid.Divide(LHSWorkspace=input_workspace, RHSWorkspace="T", OutputWorkspace=input_workspace)
    mantid.DeleteWorkspace("T")
    print "absorb done"


def generate_name_from_run(run_number, extension):
    filename = "WISH" + str(run_number)
    filename = filename + "." + extension
    print filename
    return filename


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


def shared_load_files(extension, filename, output, spectrum_max, spectrum_min, is_monitor):
    if not (extension == "nxs" or extension == "raw" or extension[0] == "s"):
        return False
    mantid.Load(Filename=filename, OutputWorkspace=output, SpectrumMin=spectrum_min,
                SpectrumMax=spectrum_max)
    mantid.Rebin(InputWorkspace=output, OutputWorkspace=output, Params='6000,-0.00063,110000')
    if not is_monitor:
        mantid.MaskBins(InputWorkspace=output,  OutputWorkspace=output, XMin=99900, XMax=106000)
    if extension == "nxs":
        mantid.ConvertUnits(InputWorkspace=output, OutputWorkspace=output, Target="Wavelength", Emode="Elastic")
    return True


def split_string(input_string):
    index_p = 0
    for i in range(len(input_string)):
        if input_string[i] == "+":
            index_p = i
    if index_p != 0:
        return int(input_string[:index_p]), int(input_string[index_p + 1:len(input_string)])


def split_string_event(input_string):
    # this assumes the form nxs_event_label_tmin_tmax
    section = input_string.split('_')
    label = section[2]
    t_min = section[3]
    t_max = section[4]
    return label, t_min, t_max


def split_workspace(focus_ws):
    for workspace_index in range(Wish.NUM_PANELS):
        out = focus_ws[:len(focus_ws) - 3] + "-" + str(workspace_index + 1) + "foc"
        mantid.ExtractSingleSpectrum(InputWorkspace=focus_ws, OutputWorkspace=out, WorkspaceIndex=workspace_index)
        mantid.DeleteWorkspace(focus_ws)
    return
