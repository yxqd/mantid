# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#

import unittest

from mantidqt.project import encoderfactory


class EncoderFactoryTest(unittest.TestCase):
    def test_can_find_encoder(self):
        encoder_factory = encoderfactory.EncoderFactory()
        try:
            encoder_factory.find_encoder("ToFConverter")
        except ValueError:
            self.fail("Value error raised by encoder factory")

    def test_throws_when_cant_find_encoder(self):
        encoder_factory = encoderfactory.EncoderFactory()
        with self.assertRaises(ValueError):
            encoder_factory.find_encoder("asihbfahufuhafuhasbfasbuasbdbasfhasfbasasfbuabfuh")
