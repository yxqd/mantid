# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#
import mantidqt.project.coders
from mantidqt.project.coders import interfaces


class EncoderFactory(object):
    def __init__(self):
        self.list_of_encoders = []

    def find_encoder(self, encoder_tag):
        """

        :param encoder_tag:
        :return:
        """
        for ii in self.list_of_encoders:
            if not isinstance(ii, interfaces.Encoder):
                raise ValueError("An invalid object was found in the list of Encoders")

            if encoder_tag in ii.get_tags():
                return ii

        raise ValueError("Invalid encoder tag passed to find_encoder")
