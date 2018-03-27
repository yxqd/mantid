#include "MantidQtWidgets/Common/Batch/JobTreeView.h"
#include <QKeyEvent>
#include <QStandardItemModel>
#include <iostream>
namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

JobTreeView::JobTreeView(QWidget *parent)
    : JobTreeView({"First", "Second", "Third", "Fourth"}, parent) {

  appendChildRowOf(RowLocation({}),
                   {"First", "Second", "ThirdItem", "FourthItem"});
  appendChildRowOf(RowLocation({}));
  insertChildRowOf(RowLocation({}), 1,
                   {"OtherItem", "THis", "THisOtherThing", "This FourthThing"});

  removeRowAt(RowLocation({}));

  // auto *root = m_model.invisibleRootItem();
  // auto *hello = new QStandardItem(QIcon(":/delete.png"), "Hello");
  // auto *world = new QStandardItem("World");
  //
  // m_model.setHorizontalHeaderLabels({"First Colunm", "Second Column"});
  // root->appendRow({hello, world});
  //
  // auto *child1 = new QStandardItem("First Child");
  // auto *child2 = new QStandardItem("Second Child");
  //
  // hello->appendRow({child1, child2});
}

QStandardItem *JobTreeView::emptyRow() const {
  return new QStandardItem(1, model()->columnCount());
}

QModelIndex JobTreeView::modelIndexAt(RowLocation const &location) const {
  auto parentIndex = QModelIndex();
  for (auto const &pathComponent : location.path())
    parentIndex = model()->index(pathComponent, 0, parentIndex);
  return parentIndex;
}

QStandardItem *JobTreeView::modelItemAt(RowLocation const &location) const {
  return modelItemFromIndex(modelIndexAt(location));
}

QStandardItem *JobTreeView::parentItemOf(RowLocation const &location) const {
  return modelItemFromIndex(parentIndexOf(location));
}

QModelIndex JobTreeView::parentIndexOf(RowLocation const &location) const {
  auto parentIndex = QModelIndex();
  auto &path = location.path();
  for (auto it = path.cbegin(); it != path.cend() - 1; ++it)
    parentIndex = model()->index(*it, 0, parentIndex);
  return parentIndex;
}

QStandardItem *JobTreeView::modelItemFromIndex(QModelIndex const &index) const {
  if (index.isValid())
    return m_model.itemFromIndex(index);
  else
    return m_model.invisibleRootItem();
}

void JobTreeView::appendChildRowOf(RowLocation const &parentLocation,
                                   std::vector<std::string> const &rowText) {
  assertOrThrow(static_cast<int>(rowText.size()) == model()->columnCount(),
                "appendChildRowOf: The number of columns in each row must be the same");
  auto *parentItem = modelItemAt(parentLocation);
  parentItem->appendRow(rowFromStringVector(rowText));
}

void JobTreeView::insertChildRowOf(RowLocation const &parentLocation,
                                   int beforeRow,
                                   std::vector<std::string> const &rowText) {
  assertOrThrow(static_cast<int>(rowText.size()) == model()->columnCount(),
                "insertChildRowOf: The number of columns in each row must be the same");
  auto *parentItem = modelItemAt(parentLocation);
  parentItem->insertRow(beforeRow, rowFromStringVector(rowText));
}

void JobTreeView::appendChildRowOf(RowLocation const &parentLocation) {
  appendChildRowOf(parentLocation,
                   std::vector<std::string>(model()->rowCount(), ""));
}

QList<QStandardItem *> JobTreeView::rowFromStringVector(
    std::vector<std::string> const &rowTextItems) const {
  auto qStringList = QList<QStandardItem *>();
  for (auto &rowText : rowTextItems)
    qStringList.append(new QStandardItem(QString::fromStdString(rowText)));
  return qStringList;
}

void JobTreeView::removeRowAt(RowLocation const &location) {
  if (location.isRoot()) {
    model()->removeRows(0, model()->rowCount());
  } else {
    auto *parent = parentItemOf(location);
    parent->removeRow(location.rowRelativeToParent());
  }
}

std::vector<std::string>
JobTreeView::rowTextAt(RowLocation const &location) const {
  return std::vector<std::string>();
}

void JobTreeView::setRowTextAt(RowLocation const &location,
                               std::vector<std::string> const &rowText) {}

std::string JobTreeView::textAt(RowLocation location, int column) {
  return std::string();
}

void JobTreeView::setTextAt(RowLocation location, int column,
                            std::string const &cellText) {}

JobTreeView::JobTreeView(QStringList const &columnHeadings, QWidget *parent)
    : QTreeView(parent), m_model(this) {
  setModel(&m_model);
  setItemDelegate(new CellDelegate(nullptr, *this));
  m_model.setHorizontalHeaderLabels(columnHeadings);
}

QModelIndex JobTreeView::addChild(QModelIndex const &parent) {
  QStandardItem *item = nullptr;
  if (parent.isValid())
    item = m_model.itemFromIndex(parent);
  else
    item = m_model.invisibleRootItem();

  auto *child1 = new QStandardItem("New");
  auto *child2 = new QStandardItem("Item");
  item->appendRow({child1, child2});

  auto lastRowIndex = item->rowCount();
  auto itemIndex = model()->index(lastRowIndex - 1, 0, parent);
  return itemIndex;
}

void JobTreeView::expand(QModelIndex const &index) {
  auto expandAt = index;
  while (expandAt.isValid()) {
    setExpanded(expandAt, true);
    expandAt = model()->parent(expandAt);
  }
}

QModelIndex JobTreeView::addSibling(QModelIndex const &index) {
  return addChild(model()->parent(index));
}

void JobTreeView::expandAndMoveTo(QModelIndex const &index) {
  clearSelection();
  expand(index);
  setCurrentIndex(index);
  edit(index);
}

QModelIndex JobTreeView::addExpandAndMoveToChild(QModelIndex const &index) {
  auto newIndex = addChild(index);
  expandAndMoveTo(newIndex);
  return newIndex;
}

void JobTreeView::keyPressEvent(QKeyEvent *event) {
  if (event->key() == Qt::Key_Return) {
    event->accept();
    if (event->modifiers() & Qt::ControlModifier)
      addExpandAndMoveToChild(currentIndex());
    else {
      auto index = currentIndex();
      auto below = QModelIndex();
      if (index.row() + 1 < model()->rowCount(index.parent()))
        below = index.sibling(index.row() + 1, 0);
      else
        below = addSibling(index);
      expandAndMoveTo(below);
    }
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
      else {
        auto newCell = addSibling(currentIndex());
        expand(newCell);
        return newCell;
      }
    }
  } else if (cursorAction == QAbstractItemView::MovePrevious) {
    auto index = currentIndex();
    if (index.isValid()) {
      if (index.column() >= 1)
        return index.sibling(index.row(), index.column() - 1);
      else if (index.row() >= 1)
        return index.sibling(index.row() - 1, model()->columnCount() - 1);
      else {
        auto parent = model()->parent(index);
        if (parent.isValid())
          return parent.sibling(parent.row(), model()->columnCount() - 1);
        else
          return QModelIndex();
      }
    }
  }
  return QTreeView::moveCursor(cursorAction, modifiers);
}

} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
