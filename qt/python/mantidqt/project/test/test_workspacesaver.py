# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#

import unittest

from os import listdir, mkdir
from os.path import isdir, expanduser
from shutil import rmtree

from mantid.api import AnalysisDataService as ADS, IMDEventWorkspace
from mantid.dataobjects import MDHistoWorkspace, MaskWorkspace
from mantidqt.project import workspacesaver
from mantid.simpleapi import CreateSampleWorkspace, CreateMDHistoWorkspace, LoadMD, LoadMask, MaskDetectors, ExtractMask
from testhelpers import WorkspaceCreationHelper

working_directory = expanduser("~") + "/workspace_saver_test_directory"


class WorkspaceSaverTest(unittest.TestCase):
    def setUp(self):
        if not isdir(working_directory):
            mkdir(working_directory)

    def tearDown(self):
        ADS.clear()
        if isdir(working_directory):
            rmtree(working_directory)

    def test_getting_workspaces_from_the_ADS(self):
        #Given
        ws1 = CreateSampleWorkspace()
        ADS.addOrReplace("ws1", ws1)
        #When
        #Then
        ws_list = workspacesaver._get_workspaces_to_save()
        self.assertEqual(ws_list[0], "ws1")
        self.assertEqual(len(ws_list), 1)

    def test_when_single_workspace_in_ADS(self):
        ws_saver = workspacesaver.WorkspaceSaver(working_directory)
        ws1 = CreateSampleWorkspace()
        ws1_name = "ws1"

        ADS.addOrReplace(ws1_name, ws1)
        ws_saver.save_workspaces()

        list_of_files = listdir(working_directory)
        self.assertEqual(len(list_of_files), 1)
        self.assertEqual(list_of_files[0], ws1_name)

    def test_when_multiple_workspaces_are_in_ADS(self):
        ws_saver = workspacesaver.WorkspaceSaver(working_directory)
        ws1 = CreateSampleWorkspace()
        ws1_name = "ws1"
        ws2 = CreateSampleWorkspace()
        ws2_name = "ws2"

        ADS.addOrReplace(ws1_name, ws1)
        ADS.addOrReplace(ws2_name, ws2)
        ws_saver.save_workspaces()

        list_of_files = listdir(working_directory)
        self.assertEqual(len(list_of_files), 2)
        self.assertEqual(list_of_files[0], ws2_name)
        self.assertEqual(list_of_files[1], ws1_name)

    def test_when_MDWorkspace_is_in_ADS(self):
        ws_saver = workspacesaver.WorkspaceSaver(working_directory)
        ws1 = CreateMDHistoWorkspace(SignalInput='1,2,3,4,5,6,7,8,9', ErrorInput='1,1,1,1,1,1,1,1,1',
                                     Dimensionality='2', Extents='-1,1,-1,1', NumberOfBins='3,3', Names='A,B',
                                     Units='U,T')
        ws1_name = "ws1"

        ADS.addOrReplace(ws1_name, ws1)
        ws_saver.save_workspaces()

        list_of_files = listdir(working_directory)
        self.assertEqual(len(list_of_files), 1)
        self.assertEqual(list_of_files[0], ws1_name)
        self._load_MDWorkspace_and_test_it(ws1_name)

    def test_when_MaskedWorkspace_is_in_ADS(self):
        ws_saver = workspacesaver.WorkspaceSaver(working_directory)
        dummy_ws = WorkspaceCreationHelper.create2DWorkspaceWithFullInstrument(2, 102, False)  # no monitors
        ADS.addOrReplace("dummy_ws", dummy_ws)
        MaskDetectors(Workspace="dummy_ws", WorkspaceIndexList=1)
        ws1_name = "ws1"
        ExtractMask(InputWorkspace="dummy_ws", OutputWorkspace=ws1_name)

        ADS.remove("dummy_ws")
        ws_saver.save_workspaces()

        list_of_files = listdir(working_directory)
        self.assertEqual(len(list_of_files), 1)
        self.assertEqual(list_of_files[0], ws1_name)
        self._load_MaskWorkspace_and_test_it(ws1_name)

    def _load_MDWorkspace_and_test_it(self, save_name):
        filename = working_directory + '/' + save_name
        ws = LoadMD(Filename=filename)
        ws_is_a_mdworkspace = isinstance(ws, IMDEventWorkspace) or isinstance(ws, MDHistoWorkspace)
        self.assertEqual(ws_is_a_mdworkspace, True)

    def _load_MaskWorkspace_and_test_it(self, save_name):
        filename = working_directory + '/' + save_name
        ws = LoadMask(InputFile=filename, Instrument="INTER")
        ws_is_a_maskworkspace = isinstance(ws, MaskWorkspace)
        self.assertEqual(ws_is_a_maskworkspace, True)
