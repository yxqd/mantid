from stresstesting import MantidStressTest
from wish.reduce import Wish

from mantid import config
import mantid.simpleapi as mantid
import os
import shutil

DIRS = config['datasearch.directories'].split(';')

# Relative to system data folder
working_folder_name = "WISH"

# Relative to working folder
input_folder_name = "input"
output_folder_name = "output"

# Relative to input folder
calibration_folder_name = "Cal"


# Generate paths for the tests
# This implies DIRS[0] is the system test data folder
working_dir = os.path.join(DIRS[0], working_folder_name)

input_dir = os.path.join(working_dir, input_folder_name)
output_dir = os.path.join(working_dir, output_folder_name)

calibration_dir = os.path.join(input_dir, calibration_folder_name)


class WISHPowderReductionTest(MantidStressTest):
    # still missing required files check with ./systemtest -R PowderReduction --showskipped
    def requiredFiles(self):
        input_files = ["WISHvana41865-1foc.nxs", "WISHvana41865-2foc.nxs", "WISHvana41865-3foc.nxs",
                       "WISHvana41865-4foc.nxs", "WISHvana41865-5foc.nxs", "WISHvana41865-6foc.nxs",
                       "WISHvana41865-7foc.nxs", "WISHvana41865-8foc.nxs", "WISHvana41865-9foc.nxs",
                       "WISHvana41865-10foc.nxs", "emptyinst38581-1foc.nxs", "emptyinst38581-2foc.nxs",
                       "emptyinst38581-3foc.nxs", "emptyinst38581-4foc.nxs", "emptyinst38581-5foc.nxs",
                       "emptyinst38581-6foc.nxs", "emptyinst38581-7foc.nxs", "emptyinst38581-8foc.nxs",
                       "emptyinst38581-9foc.nxs", "emptyinst38581-10foc.nxs"]

        input_files = [os.path.join(calibration_dir, files) for files in input_files]
        return input_files

    def cleanup(self):
        shutil.rmtree(output_dir)

    def runTest(self):
        os.makedirs(output_dir)
        wish_test = Wish("__main__", calibration_dir+"/", input_dir, output_dir, True)
        wish_test.run_script(41870)
        self.clearWorkspaces()

    def validate(self):
        return "w41870_2-9foc", "WISH41870_2-9raw.nxs", \
               "w41870_3-8foc", "WISH41870_3-8raw.nxs", \
               "w41870_4-7foc", "WISH41870_4-7raw.nxs", \
               "w41870_5-6foc", "WISH41870_5-6raw.nxs"

    def clearWorkspaces(self):
        deletews = ["w41870_" + str(i) + "foc" for i in range(1, 11)]
        for ws in deletews:
            mantid.DeleteWorkspace(ws)
            mantid.DeleteWorkspace(ws + "-d")

    # Skip test when on builds as extremely slow, run only as reversion test for wish script
    def skipTests(self):
        return True