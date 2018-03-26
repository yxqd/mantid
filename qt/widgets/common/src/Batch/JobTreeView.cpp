#include "MantidQtWidgets/Common/Batch/JobTreeView.h"
namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

JobTreeView::JobTreeView(QWidget *parent)
    : QTreeView(parent), m_modelAdapter(this) {
  setModel(&m_modelAdapter);
  setItemDelegate(new CellDelegate(nullptr, m_modelAdapter));
}

QModelIndex JobTreeView::moveCursor(CursorAction cursorAction,
                                    Qt::KeyboardModifiers modifiers) {
  if (cursorAction == QAbstractItemView::MoveNext) {
    auto index = currentIndex();
    if (index.isValid()) {
      if (index.column() + 1 < model()->columnCount())
        return index.sibling(index.row(), index.column() + 1);
      else if (index.row() + 1 < model()->rowCount(index.parent()))
        return index.sibling(index.row() + 1, 0);
      else
        return QModelIndex();
    }
  } else if (cursorAction == QAbstractItemView::MovePrevious) {
    auto index = currentIndex();
    if (index.isValid()) {
      if (index.column() >= 1)
        return index.sibling(index.row(), index.column() - 1);
      else if (index.row() >= 1)
        return index.sibling(index.row() - 1, model()->columnCount() - 1);
      else
        return QModelIndex();
    }
  }

  return QTreeView::moveCursor(cursorAction, modifiers);
}

} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
