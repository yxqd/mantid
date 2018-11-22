# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#
from json import dumps
from qtpy.QtWidgets import QFileDialog, QApplication
from os import makedirs
from os.path import isdir

from mantidqt.project import workspacesaver, encoderfactory
from mantidqt.io import open_a_file_dialog
from mantid import logger


# Static method for getting the save location
def _get_save_location():
    return open_a_file_dialog(accept_mode=QFileDialog.AcceptSave, file_mode=QFileDialog.DirectoryOnly)


ENCODED_FILE_NAME = "mantidsave.interfaces"


class ProjectSaver(object):
    def __init__(self, save_location):
        """

        :param save_location: String; the place to save your project to
        """
        self.encoderfactory = encoderfactory.EncoderFactory()
        self.directory = save_location

    def save_project(self, workspace_to_save=None, interfaces_to_save=None):
        """

        :param workspace_to_save: List; of Strings that will have workspace names in it, if None will save all
        :param interfaces_to_save: List; of Strings that will have interface tags in it, if None will save all
        :return:
        """
        # Check this isn't saving a blank project file
        if workspace_to_save is None and interfaces_to_save is None:
            logger.warning("Can not save an empty project")
            return

        # Save workspaces to that location
        workspace_saver = workspacesaver.WorkspaceSaver(directory=self.directory,
                                                        workspaces_to_save=workspace_to_save)
        workspace_saver.save_workspaces()

        # Get interface details in dicts
        dictionaries_to_save = self._encode_interfaces(interfaces_to_save)

        # Pass dicts to Project Writer
        writer = ProjectWriter(dictionaries_to_save, self.directory, workspace_saver.get_output_list())
        writer.write_out()

    def _encode_interfaces(self, interfaces_to_save=None):
        encoders = self._get_encoders(interfaces_to_save)
        if encoders is None:
            return {}
        dicts = {}
        for ii in encoders:
            # Save using the first tag of the interface
            dicts[ii.get_tags()[0]] = ii.encode()
        return dicts

    def _get_encoders(self, interfaces_to_save=None):
        list_of_encoders = []
        if interfaces_to_save is None:
            return None
        for ii in interfaces_to_save:
            list_of_encoders.append(self.encoderfactory.find_encoder(ii))
        return list_of_encoders


# Static private method to create JSON from objects and return a string
def _create_out_string_json(dictionary):
    return dumps(dictionary)


class ProjectWriter(object):
    def __init__(self, dicts, save_location, workspace_names):
        self.dicts_to_save = dicts
        self.workspace_names = workspace_names
        self.directory = save_location

    def write_out(self):
        # Get the JSON string versions
        workspace_interface_dict = {"Workspaces": self.workspace_names, "Interfaces": self.dicts_to_save}
        json_string = _create_out_string_json(workspace_interface_dict)

        # Open file and save the string to it alongside the workspace_names
        file_name = self.directory + '/' + ENCODED_FILE_NAME
        if not isdir(self.directory):
            makedirs(self.directory)
        f = open(file_name, 'w+')
        f.write(json_string)
