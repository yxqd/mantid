#ifndef MANTIDQTCUSTOMINTERFACES_IDATACOMPARISONPRESENTER_H_
#define MANTIDQTCUSTOMINTERFACES_IDATACOMPARISONPRESENTER_H_

namespace MantidQt {
namespace CustomInterfaces {

/**
Interface for the presenter as in the MVP pattern. Methods defined here will be
used by the view to notify the presenter that something happened.

Copyright &copy; 2016 ISIS Rutherford Appleton Laboratory, NScD
Oak Ridge National Laboratory & European Spallation Source

This file is part of Mantid.

Mantid is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Mantid is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

File change history is stored at: <https://github.com/mantidproject/mantid>
Code Documentation is available at: <http://doxygen.mantidproject.org>
*/

class IDataComparisonPresenter {
public:
  virtual ~IDataComparisonPresenter(){};

  /// User actions triggered from the (passive) view that need handling by the
  /// presenter
  enum Notification {
    AddWorkspace,
    PlotWorkspaces,
    PlotDiffWorkspaces,
    RemoveAllWorkspaces,
    RemoveSelectedWorkspaces,
    RemoveDiffWorkspace,
	WorkspaceIndexChanged,
	ColorChanged
  };

  virtual void notify(IDataComparisonPresenter::Notification notification) = 0;
};
}
}

#endif /* MANTIDQTCUSTOMINTERFACES_IDATACOMPARISONPRESENTER_H_ */
