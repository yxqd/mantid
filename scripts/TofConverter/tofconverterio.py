# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#

from mantidqt.project.coders import interfaces
from mantid.kernel import logger


INTERFACE_KEYS = ["TofConverter"]


class ToFConverterEncoder(interfaces.Encoder):
    def encode(self, obj):
        return {}

    def get_interface_keys(self):
        return INTERFACE_KEYS


class ToFConverterDecoder(interfaces.Decoder):
    def decode(self, obj):
        if obj is not {}:
            logger.warning("Data was loaded for the TofConverter that should not exist, check the encoder")

        return

    def get_interface_keys(self):
        return INTERFACE_KEYS
