#ifndef MANTIDQTCUSTOMINTERFACES_DATACOMPARISONPRESENTER_H_
#define MANTIDQTCUSTOMINTERFACES_DATACOMPARISONPRESENTER_H_

#include "MantidQtCustomInterfaces/DataComparison/IDataComparisonPresenter.h"
#include "MantidQtAPI/WorkspaceObserver.h"

namespace MantidQt {
namespace CustomInterfaces {

class IDataComparisonView;
/**
Concrete implementation of the presenter.

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

class DataComparisonPresenter : public IDataComparisonPresenter,
                                public MantidQt::API::WorkspaceObserver {

public:
  /// Constructor
  DataComparisonPresenter(IDataComparisonView *view);
  /// Destructor
  ~DataComparisonPresenter() override;
  /// Notifications
  void notify(IDataComparisonPresenter::Notification notification) override;

private:
  /// The view
  IDataComparisonView *m_view;
  /// The diff workspace
  std::string m_diffWsName;

  /// Add workspace to table
  void addWorkspace(Mantid::API::Workspace_const_sptr ws);
  /// Add workspace to interface
  void addWorkspace();
  /// Plot workspaces
  void plotWorkspaces();
  /// Plot diff workspace
  void plotDiffWorkspace();
  /// Remove all workspaces
  void removeAllWorkspaces();
  /// Remove selected workspaces
  void removeSelectedWorkspaces();
  /// Remove diff workspace
  void removeDiffWorkspace();
  /// Get initial color for new ws
  int getInitialColourIndex();

  // Handlers for ADS events
  void
  preDeleteHandle(const std::string &wsName,
                  const boost::shared_ptr<Mantid::API::Workspace> ws) override;
  void renameHandle(const std::string &oldName,
                    const std::string &newName) override;
  void afterReplaceHandle(
      const std::string &wsName,
      const boost::shared_ptr<Mantid::API::Workspace> ws) override;
};
}
}
#endif /* MANTIDQTCUSTOMINTERFACES_DATACOMPARISONPRESENTER_H_ */
