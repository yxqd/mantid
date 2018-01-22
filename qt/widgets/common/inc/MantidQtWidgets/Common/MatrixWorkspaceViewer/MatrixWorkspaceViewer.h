#ifndef MANTID_WIDGETS_MATRIX_WORKSPACE_VIEWER_H
#define MANTID_WIDGETS_MATRIX_WORKSPACE_VIEWER_H

#include "MantidMatrixModel.h"

#include "MantidQtWidgets/Common/DllOption.h"
#include "MantidQtWidgets/Common/WorkspaceObserver.h"

#include "MantidAPI/MatrixWorkspace_fwd.h"

#include <QFrame>
#include <QPixmap>
#include <QPointer>
#include <QTabWidget>
#include <QTableView>

namespace MantidQt {
namespace MantidWidgets {

/** MantidMatrix is the class that represents a Qtiplot window for displaying
workspaces.
It has separate tabs for displaying spectrum values, bin boundaries, and errors.

@author Roman Tolchenov, Tessella Support Services plc

Copyright &copy; 2007 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
National Laboratory & European Spallation Source

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

File change history is stored at: <https://github.com/mantidproject/mantid>.
Code Documentation is available at: <http://doxygen.mantidproject.org>

*/
class EXPORT_OPT_MANTIDQT_COMMON MatrixWorkspaceViewer : public QFrame, MantidQt::API::WorkspaceObserver {
  Q_OBJECT
  MatrixWorkspaceViewer(const MatrixWorkspaceViewer&) = delete;

public:
  MatrixWorkspaceViewer(QWidget *parent=nullptr);
  void setWorkspace(Mantid::API::MatrixWorkspace_sptr ws, int start = -1, int end = -1);
  void setName(QString name) {m_name = name;}
  int numRows() const { return m_rows; }
  int numCols() const { return m_cols; }

private:
  void setup(Mantid::API::MatrixWorkspace_sptr ws, int start, int end);
  void connectTableView(QTableView *, MantidMatrixModel *);

  QString m_label;
  QString m_name;
  Mantid::API::MatrixWorkspace_sptr m_workspace;
  QTabWidget *m_tabs;
  QTableView *m_table_viewY;
  QTableView *m_table_viewX;
  QTableView *m_table_viewE;
  QPointer<MantidMatrixModel> m_modelY;
  QPointer<MantidMatrixModel> m_modelX;
  QPointer<MantidMatrixModel> m_modelE;
  QColor m_bk_color;
  QPixmap m_matrix_icon;
  // The tab labels
  QString m_YTabLabel, m_XTabLabel, m_ETabLabel;
  int m_column_width;
  // index to identify the previous view on tab switch
  int m_PrevIndex;

  double x_start, //!< X value corresponding to column 1
      x_end,      //!< X value corresponding to the last column
      y_start,    //!< Y value corresponding to row 1
      y_end;      //!< Y value corresponding to the last row
  int m_rows, m_cols;
  int m_startRow;
  int m_endRow;
  int m_workspaceTotalHist;
  bool m_histogram;
  double m_min;           // Saved minimum Y-value
  double m_max;           // Saved maximum Y-value
  bool m_are_min_max_set; // If true ::range does not iterate over WS to find
                          // min and max but uses m_min and m_max instead


};

} // namespace MantidQt
} // namespace MantidWidgets

#endif // MANTID_WIDGETS_MATRIX_WORKSPACE_VIEWER_H
