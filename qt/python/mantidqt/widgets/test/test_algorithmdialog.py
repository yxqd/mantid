from mock import Mock

from mantid.api import AlgorithmFactory, AlgorithmManager, PythonAlgorithm
from mantid.kernel import Direction, FloatArrayProperty
from mantidqt.dialogs.genericdialog import GenericDialog
from mantidqt.utils.qt.test.gui_window_test import GuiWindowTest

mock_alg_executed = Mock()
MOCK_PARAMETER = 42.42


class AlgorithmDialogMockAlgorithm(PythonAlgorithm):

    def category(self):
        return 'Examples'

    def PyInit(self):
        self.declareProperty("InValue", 0)
        self.declareProperty("DoubleValue", 1.0)
        self.declareProperty(FloatArrayProperty("Floats", direction=Direction.Input))
        self.declareProperty("OutValue", 0, direction=Direction.Output)

    def PyExec(self):
        self.debug()
        global mock_alg_executed
        mock_alg_executed(MOCK_PARAMETER)


class TestAlgorithmDialogVisual(GuiWindowTest):
    def setUp(self):
        AlgorithmFactory.subscribe(AlgorithmDialogMockAlgorithm)
        global mock_alg_executed
        mock_alg_executed = Mock()

    def tearDown(self):
        AlgorithmFactory.unsubscribe('AlgorithmDialogMockAlgorithm', 1)
        global mock_alg_executed
        mock_alg_executed = None

    def create_widget(self):
        dialog = GenericDialog()
        alg = AlgorithmManager.create('AlgorithmDialogMockAlgorithm')
        dialog.setAlgorithm(alg)
        dialog.initializeLayout()
        dialog.show()
        return dialog

    def test_closes(self):
        self.click_button_by_text("Close")
        self.assertEqual(0, mock_alg_executed.call_count)

    def test_runs(self):
        self.click_button_by_text("Run")
        global mock_alg_executed
        mock_alg_executed.assert_called_once_with(MOCK_PARAMETER)

    # def test_fields_shown(self):
    #     self.get_child_by_text(QLabel, "InValue")
    #     self.get_child_by_text(QLabel, "DoubleValue")
    #     self.get_child_by_text(QLabel, "Floats")

def debug():
    PYTHON_ROOT = "c:/users/qbr77747/apps/miniconda3"
    import os
    import sys
    sys.path.append(os.path.join(PYTHON_ROOT, "lib/site-packages"))
    import pydevd
    pydevd.settrace('localhost', port=44444, stdoutToServer=True, stderrToServer=True)
