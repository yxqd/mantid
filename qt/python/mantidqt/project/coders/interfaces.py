# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#


class Encoder(object):
    """
    This is an abstract class to be implemented to allow integration with the encoder factory
    """
    def encode(self, obj):
        """
        This is the method that should encode the interface
        :return: Return a dictionary of information relevant to the class
        """
        raise NotImplementedError("Class %s doesn't implement encode()" % self.__class__.__name__)

    def get_interface_keys(self):
        """

        :return:
        """
        raise NotImplementedError("Class %s doesn't implement get_tags()" % self.__class__.__name__)


class Decoder(object):
    """
    This is an abstract class to be implemented to allow integration with the decoder factory
    """
    def decode(self, obj):
        raise NotImplementedError("Class %s doesn't implement decode()" % self.__class__.__name__)

    def get_interface_keys(self):
        raise NotImplementedError("Class %s doesn't implement get_tags()" % self.__class__.__name__)


class InterfaceAttributes(object):
    """
    This is an abstract class to be implemented to allow common data to be shared between the decoder and the encoder
    """
    def get_interface_keys(self):
        """

        :return:
        """
        raise NotImplementedError("Class %s doesn't implement decode()" % self.__class__.__name__)
