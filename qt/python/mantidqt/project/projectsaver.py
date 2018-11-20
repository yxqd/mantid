# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#
from json import dumps
from qtpy.QtWidgets import QFileDialog

from mantidqt.project import workspacesaver, encoderfactory
from mantidqt.io import open_a_file_dialog


# Static method for getting the save location
def _get_save_location():
    return open_a_file_dialog(accept_mode=QFileDialog.AcceptSave, file_mode=QFileDialog.DirectoryOnly)


class ProjectSaver(object):
    def __init__(self, save_location, workspace_to_save=None, interfaces_to_save=None):
        self.encoderfactory = encoderfactory.EncoderFactory()
        self.directory = save_location
        self.workspace_to_save = workspace_to_save
        self.interfaces_to_save = interfaces_to_save

    def save_project(self):
        # Save workspaces to that location
        workspace_saver = workspacesaver.WorkspaceSaver(directory=self.directory,
                                                        workspaces_to_save=self.workspace_to_save)
        workspace_saver.save_workspaces()

        # Get interface details in dicts
        dictionaries_to_save = self._encode_interfaces()

        # Pass dicts to Project Writer
        writer = ProjectWriter(dictionaries_to_save, self.directory, workspace_saver.get_output_list())
        writer.write_out()

    def _encode_interfaces(self):
        encoders = self._get_encoders()
        dicts = []
        for ii in encoders:
            dicts.append(ii.encode())
        return dicts

    def _get_encoders(self):
        list_of_encoders = []
        for ii in self.interfaces_to_save:
            list_of_encoders.append(self.encoderfactory.find_encoder(ii))
        return list_of_encoders


class ProjectWriter(object):
    def __init__(self, dicts, save_location, workspace_names):
        self.dicts_to_save = dicts
        self.workspace_names = workspace_names
        self.directory = save_location

    def write_out(self):
        # Get the JSON string versions
        json_string = self._create_out_string()

        # Open file and save the string to it alongside the workspace_names

    def _create_out_string(self):
        return dumps(self.dicts_to_save)
