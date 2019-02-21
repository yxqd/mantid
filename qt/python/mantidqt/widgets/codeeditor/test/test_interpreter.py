# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2019 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#    This file is part of the mantid workbench.
#
#
from __future__ import (absolute_import, unicode_literals)

import unittest

import six

from mantidqt.utils.qt.test import GuiTest
from mantidqt.widgets.codeeditor.interpreter import PythonFileInterpreter

if six.PY2:
    import mock
else:
    from unittest import mock


class PythonFileInterpreterTest(GuiTest):

    def test_construction(self):
        w = PythonFileInterpreter(None)
        self.assertTrue("Status: Idle", w.status.currentMessage())

    def test_empty_code_does_nothing_on_exec(self):
        w = PythonFileInterpreter(None)
        w._presenter.model.execute_async = mock.MagicMock()
        w.execute_async()
        w._presenter.model.execute_async.assert_not_called()
        self.assertTrue("Status: Idle", w.status.currentMessage())

    def test_constructor_populates_editor_with_content(self):
        w = PythonFileInterpreter(None, content='# my funky code')
        self.assertEqual('# my funky code', w.editor.text())

    def test_constructor_respects_filename(self):
        w = PythonFileInterpreter(None, filename='test.py')
        self.assertEqual('test.py', w.filename)

    def test_successful_execution(self):
        w = PythonFileInterpreter(None)
        w.editor.setText("x = 1 + 2")
        w.execute_async()
        self.assertTrue("Status: Idle", w.status.currentMessage())

    def test_clear_key_binding(self):
        test_cases = {'Ctrl+A': None, 'Shift+A': ValueError,
                      'Ctrl+AAA': ValueError, 'Ctrl+Shift+A': ValueError}
        w = PythonFileInterpreter(None)
        for key_combo, expected_result in test_cases.items():
            fail_msg = ("Failed on case '{}' with expected result '{}'"
                        "".format(key_combo, expected_result))
            if expected_result is ValueError:
                with self.assertRaises(expected_result, msg=fail_msg):
                    w.clear_key_binding(key_combo)
            else:
                self.assertEqual(w.clear_key_binding(key_combo), None,
                                 msg=fail_msg)


if __name__ == '__main__':
    unittest.main()
