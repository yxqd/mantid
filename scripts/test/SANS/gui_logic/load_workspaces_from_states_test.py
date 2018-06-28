from __future__ import (absolute_import, division, print_function)

import unittest
import sys

from mantid.kernel import config
from mantid.kernel import PropertyManagerDataService

from sans.gui_logic.presenter.run_tab_presenter import RunTabPresenter
from sans.common.enums import (SANSFacility, ReductionDimensionality, SaveType, ISISReductionMode,
                               RangeStepType, FitType)
from sans.test_helper.user_file_test_helper import (create_user_file, sample_user_file, sample_user_file_gravity_OFF)
from sans.test_helper.mock_objects import (create_mock_view)
from sans.test_helper.common import (remove_file)
from sans.common.enums import BatchReductionEntry
from sans.gui_logic.models.load_workspaces_from_states import load_workspaces_from_states


if sys.version_info.major == 3:
    from unittest import mock
else:
    import mock

class LoadWorkspacesFromStatesTest(unittest.TestCase):
    @mock.patch('sans.gui_logic.models.load_workspaces_from_states.load_workspaces')
    def test_calls_on_load_finished_after_succesful_execution(self, load_workspaces_mock):
        load_workspaces_mock.return_value = True
        presenter = mock.MagicMock()
        states = mock.MagicMock()

        load_workspaces_from_states(states, presenter)

        presenter.on_load_finished.assert_called_once_with(True)

    @mock.patch('sans.gui_logic.models.load_workspaces_from_states.load_workspaces')
    def test_calls_on_load_error_after_unsuccesful_execution(self, load_workspaces_mock):
        load_workspaces_mock.side_effect = RuntimeError('There was an error')
        load_workspaces_mock.return_value = True
        presenter = mock.MagicMock()
        states = mock.MagicMock()

        load_workspaces_from_states(states, presenter)

        presenter.on_load_error.assert_called_once_with('There was an error')


if __name__ == '__main__':
    unittest.main()