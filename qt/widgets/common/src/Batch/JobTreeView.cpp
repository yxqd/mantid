#include "MantidQtWidgets/Common/Batch/JobTreeView.h"
#include <QKeyEvent>
#include <iostream>
namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

JobTreeView::JobTreeView(QWidget *parent)
    : QTreeView(parent), m_modelAdapter(this) {
  setModel(&m_modelAdapter);
  setItemDelegate(new CellDelegate(nullptr, m_modelAdapter));
  setAlternatingRowColors(false);
}

void JobTreeView::keyPressEvent(QKeyEvent *event) {
  if (event->key() == Qt::Key_Return &&
      (event->modifiers() & Qt::ControlModifier)) {
    event->accept();
    auto parent = m_modelAdapter.parent(selectedIndexes()[0]);
    setCurrentIndex(parent);
    edit(parent);
    //auto index = m_modelAdapter.insertRow(
        //selectedIndexes()[0], QtTreeRow({QtTreeCell("New"), QtTreeCell("Row")}));
    //
    //    auto expandAt = index;
    //    while (expandAt.isValid()) {
    //      setExpanded(expandAt, true);
    //      expandAt = model()->parent(expandAt);
    //    }
    //    setCurrentIndex(index);
    //    edit(index);
  } else {
    QTreeView::keyPressEvent(event);
  }
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
