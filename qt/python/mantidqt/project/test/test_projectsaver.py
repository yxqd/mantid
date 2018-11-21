# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#
import unittest

from os import listdir
from os.path import isdir, expanduser
from shutil import rmtree

from mantid.api import AnalysisDataService as ADS
from mantid.simpleapi import CreateSampleWorkspace
from mantidqt.project import projectsaver


working_directory = expanduser("~") + "/project_saver_test"

class ProjectSaverTest(unittest.TestCase):
    def tearDown(self):
        ADS.clear()
        if isdir(working_directory):
            rmtree(working_directory)

    def test_only_one_workspace_saving(self):
        ws1_name = "ws1"
        ADS.addOrReplace(ws1_name, CreateSampleWorkspace(OutputWorkspace=ws1_name))
        project_saver = projectsaver.ProjectSaver(save_location=working_directory, save_all=True)
        file_name = working_directory + "/" + projectsaver.ENCODED_FILE_NAME
        saved_file = "[]\n[\"ws1\"]"

        project_saver.save_project()

        # Check project file is saved correctly
        f = open(file_name, "r")
        self.assertEqual(f.read(), saved_file)

        # Check workspace is saved
        list_of_files = listdir(working_directory)
        self.assertEqual(len(list_of_files), 2)
        self.assertEqual(list_of_files[0], projectsaver.ENCODED_FILE_NAME)
        self.assertEqual(list_of_files[1], ws1_name)

    #def test_only_multiple_workspaces_saving(self):

    #def test_only_saving_one_workspace_when_multiple_are_present_in_the_ADS(self):

    #def test_only_one_interface_saving(self):

    #def test_only_multiple_interfaces_saving(self):

    #def test_one_workspace_and_one_interface_saving(self):

    #def test_multiple_workspaces_and_multiple_interfaces(self):

    #def test_get_encoders_retrieves_correct_encoder(self):

    #def test_encode_interfaces_on_one_interface(self):

    #def test_encode_interfaces_on_multiple_interfaces(self):


class ProjectWriterTest(unittest.TestCase):
    def tearDown(self):
        ADS.clear()
        if isdir(working_directory):
            rmtree(working_directory)

    def test_create_out_string_on_small_dict(self):
        small_dict = {"interface1" :{"value1" : 2, "value2" : 3}, "interface2" : {"value3" : 4, "value4" : 5}}
        json_dict = "{\"interface1\": {\"value2\": 3, \"value1\": 2}, \"interface2\": {\"value4\": 5, \"value3\": 4}}"
        project_writer = projectsaver.ProjectWriter(small_dict, working_directory, [])

        output = project_writer._create_out_string()
        self.assertEqual(output, json_dict)

    def test_create_out_string_on_larger_dict(self):
        large_dict = {"interface1": {"value1": 2, "value2": 3}, "interface2": {"value3": 4, "value4": 5},
                      "interface3": {"value1": 2, "value2": 3}, "interface4": {"value3": 4, "value4": 5},
                      "interface5": {"value1": 2, "value2": 3}, "interface6": {"value3": 4, "value4": 5}}
        json_dict = "{\"interface1\": {\"value2\": 3, \"value1\": 2}, \"interface3\": {\"value2\": 3, \"value1\": 2}, "\
                    "\"interface2\": {\"value4\": 5, \"value3\": 4}, \"interface5\": {\"value2\": 3, \"value1\": 2}, "\
                    "\"interface4\": {\"value4\": 5, \"value3\": 4}, \"interface6\": {\"value4\": 5, \"value3\": 4}}"
        project_writer = projectsaver.ProjectWriter(large_dict, working_directory, [])

        output = project_writer._create_out_string()
        self.assertEqual(output, json_dict)

    def test_create_out_string_on_empty_dict(self):
        empty_dict = {}
        json_dict = "{}"
        project_writer = projectsaver.ProjectWriter(empty_dict, working_directory, [])

        output = project_writer._create_out_string()
        self.assertEqual(output, json_dict)

    def test_create_workspace_json_string(self):
        workspace_list = {"workspace_names" : ["ws1", "ws2", "ws3", "ws4"]}
        json_list = "{\"workspace_names\": [\"ws1\", \"ws2\", \"ws3\", \"ws4\"]}"
        project_writer = projectsaver.ProjectWriter({}, working_directory, workspace_list)

        output = project_writer._create_workspace_json_string()
        self.assertEqual(output, json_list)

    def test_create_workspace_json_string_with_no_workspaces(self):
        workspace_list = {"workspace_names": []}
        json_list = "{\"workspace_names\": []}"
        project_writer = projectsaver.ProjectWriter({}, working_directory, workspace_list)

        output = project_writer._create_workspace_json_string()
        self.assertEqual(output, json_list)

    def test_write_out_on_just_dicts(self):
        workspace_list = {"workspace_names": []}
        small_dict = {"interface1": {"value1": 2, "value2": 3}, "interface2": {"value3": 4, "value4": 5}}
        project_writer = projectsaver.ProjectWriter(small_dict, working_directory, workspace_list)
        file_name = working_directory + "/" + projectsaver.ENCODED_FILE_NAME
        saved_file = "{\"interface1\": {\"value2\": 3, \"value1\": 2}, \"interface2\": {\"value4\": 5, \"value3\": 4}" \
                     "}\n{\"workspace_names\": []}"

        project_writer.write_out()

        f = open(file_name, "r")
        self.assertEqual(f.read(), saved_file)

    def test_write_out_on_just_workspaces(self):
        workspace_list = {"workspace_names": ["ws1", "ws2", "ws3", "ws4"]}
        small_dict = {}
        project_writer = projectsaver.ProjectWriter(small_dict, working_directory, workspace_list)
        file_name = working_directory + "/" + projectsaver.ENCODED_FILE_NAME
        saved_file = "{}\n{\"workspace_names\": [\"ws1\", \"ws2\", \"ws3\", \"ws4\"]}"

        project_writer.write_out()

        f = open(file_name, "r")
        self.assertEqual(f.read(), saved_file)

    def test_write_out_on_both_workspaces_and_dicts(self):
        workspace_list = {"workspace_names": ["ws1", "ws2", "ws3", "ws4"]}
        small_dict = {"interface1": {"value1": 2, "value2": 3}, "interface2": {"value3": 4, "value4": 5}}
        project_writer = projectsaver.ProjectWriter(small_dict, working_directory, workspace_list)
        file_name = working_directory + "/" + projectsaver.ENCODED_FILE_NAME
        saved_file = "{\"interface1\": {\"value2\": 3, \"value1\": 2}, \"interface2\": {\"value4\": 5, \"value3\": 4}" \
                     "}\n{\"workspace_names\": [\"ws1\", \"ws2\", \"ws3\", \"ws4\"]}"

        project_writer.write_out()

        f = open(file_name, "r")
        self.assertEqual(f.read(), saved_file)