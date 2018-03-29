#include "MantidQtWidgets/Common/Batch/JobTreeView.h"
#include "MantidQtWidgets/Common/Batch/QtTreeCursorNavigation.h"
#include <QKeyEvent>
#include <QStandardItemModel>
#include <QSortFilterProxyModel>
#include <iostream>
namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

class QtFilterLeafNodes : public QSortFilterProxyModel {
public:
  QtFilterLeafNodes(QObject *parent = nullptr) : QSortFilterProxyModel(parent){};

protected:
  bool filterAcceptsRow(int row, const QModelIndex &parent) const override {
    QModelIndex index = sourceModel()->index(row, 0, parent);

    if (index.isValid()) {
      if (index.data().toString().contains("App"))
        return true;

      int rows = sourceModel()->rowCount(index);
      for (auto r = 0; r < rows; r++)
        if (filterAcceptsRow(r, index))
          return true;
      return false;
    } else {
      return false;
    }
  }
};

JobTreeView::JobTreeView(QWidget *parent)
    : JobTreeView({"Runs", "Angle", "Transmission Runs", "Q min", "Q max", "dQ/Q", "Scale", "Options"}, parent) {

  for(auto i = 0; i < model()->columnCount(); ++i)
    resizeColumnToContents(i);

  appendChildRowOf(RowLocation({}));

  

}

QModelIndex JobTreeView::modelIndexAt(RowLocation const &location,
                                      int column) const {
  auto parentIndex = QModelIndex();
  if (!location.isRoot()) {
    auto &pathComponents = location.path();
    auto it = pathComponents.cbegin();
    for (; it != pathComponents.cend() - 1; ++it)
      parentIndex = model()->index(*it, column, parentIndex);
    assertOrThrow(
        model()->hasIndex(location.rowRelativeToParent(), 0, parentIndex),
        "modelItemAt: Location refers to an index which does not "
        "exist in the model.");
    return model()->index(location.rowRelativeToParent(), column, parentIndex);
  } else
    return parentIndex;
}

void JobTreeView::subscribe(IJobTreeViewSubscriber &subscriber) {
  m_subscriber = &subscriber;
}

RowLocation JobTreeView::rowLocationAt(QModelIndex const &index) const {
  if (index.isValid()) {
    auto pathComponents = RowPath();
    auto currentIndex = index;
    while(index.isValid()) {
      pathComponents.insert(pathComponents.begin(), currentIndex.row());
      currentIndex = index.parent();
    }
    return RowLocation(pathComponents);
  } else {
    return RowLocation({});
  }
}

QStandardItem *JobTreeView::modelItemAt(RowLocation const &location,
                                        int column) const {
  return modelItemFromIndex(modelIndexAt(location, column));
}

QStandardItem *JobTreeView::parentItemOf(RowLocation const &location,
                                         int column) const {
  return modelItemFromIndex(parentIndexOf(location, column));
}

QModelIndex JobTreeView::parentIndexOf(RowLocation const &location,
                                       int column) const {
  auto parentIndex = QModelIndex();
  auto &path = location.path();
  for (auto it = path.cbegin(); it != path.cend() - 1; ++it)
    parentIndex = model()->index(*it, column, parentIndex);
  return parentIndex;
}

QStandardItem *JobTreeView::modelItemFromIndex(QModelIndex const &index) const {
  if (index.isValid()) {
    auto *item = m_model.itemFromIndex(index);
    assertOrThrow(item != nullptr,
                  "modelItemFromIndex: Index must point to a valid item.");
    return item;
  } else
    return m_model.invisibleRootItem();
}

void JobTreeView::appendChildRowOf(RowLocation const &parentLocation,
                                   std::vector<std::string> const &rowText) {
  assertOrThrow(
      static_cast<int>(rowText.size()) == model()->columnCount(),
      "appendChildRowOf: The number of columns in each row must be the same");
  make(appendedChildRow(modelIndexAt(parentLocation), rowFromRowText(rowText)));
}

void JobTreeView::insertChildRowOf(RowLocation const &parentLocation,
                                   int beforeRow,
                                   std::vector<std::string> const &rowText) {
  assertOrThrow(
      static_cast<int>(rowText.size()) == model()->columnCount(),
      "insertChildRowOf: The number of columns in each row must be the same");
  make(insertedChildRow(modelIndexAt(parentLocation), beforeRow,
                        rowFromRowText(rowText)));
}

void JobTreeView::appendChildRowOf(RowLocation const &parentLocation) {
  make(appendedEmptyChildRow(modelIndexAt(parentLocation)));
}

QList<QStandardItem *>
JobTreeView::rowFromRowText(std::vector<std::string> const &rowText) const {
  auto rowCells = QList<QStandardItem *>();
  for (auto &cellText : rowText)
    rowCells.append(new QStandardItem(QString::fromStdString(cellText)));
  return rowCells;
}

std::vector<std::string>
JobTreeView::rowTextFromRow(QModelIndex firstCellIndex) const {
  auto rowText = std::vector<std::string>();
  rowText.reserve(model()->columnCount());

  for (auto i = 0; i < model()->columnCount(); i++) {
    auto *cell =
        modelItemFromIndex(firstCellIndex.sibling(firstCellIndex.row(), i));
    rowText.emplace_back(cell->text().toStdString());
  }
  return rowText;
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
  return rowTextFromRow(modelIndexAt(location));
}

void JobTreeView::setRowTextAt(RowLocation const &location,
                               std::vector<std::string> const &rowText) {
  assertOrThrow(static_cast<int>(rowText.size()) == model()->columnCount(),
                "Replacement row text must have the same number of columns as "
                "the original.");
  assertOrThrow(!location.isRoot(), "Cannot set the text for the root node. "
                                    "Did you mean to use setHeaderText "
                                    "instead?");
  auto firstCellIndex = modelIndexAt(location);
  for (auto i = 0; i < model()->columnCount(); i++) {
    auto *cell = modelItemFromIndex(
        firstCellIndex.sibling(location.rowRelativeToParent(), i));
    auto replacementText = QString::fromStdString(rowText[i]);
    cell->setText(replacementText);
  }
}

std::string JobTreeView::textAt(RowLocation location, int column) {
  auto *cell = modelItemAt(location, column);
  return cell->text().toStdString();
}

void JobTreeView::setTextAt(RowLocation location, int column,
                            std::string const &cellText) {
  auto *cell = modelItemAt(location, column);
  cell->setText(QString::fromStdString(cellText));
}

JobTreeView::JobTreeView(QStringList const &columnHeadings, QWidget *parent)
    : QTreeView(parent), m_model(this) {
  setModel(&m_model);
  setItemDelegate(new CellDelegate(nullptr, *this));
  m_model.setHorizontalHeaderLabels(columnHeadings);
}

QModelIndex JobTreeView::appendedChildRow(QModelIndex const &parent,
                                          QList<QStandardItem *> cells) {
  auto *const parentItem = modelItemFromIndex(parent);
  parentItem->appendRow(cells);
  return navigation().lastRowInThisNode(parent);
}

QModelIndex JobTreeView::appendedEmptyChildRow(QModelIndex const &parent) {
  auto *const parentItem = modelItemFromIndex(parent);

  auto cells = QList<QStandardItem *>();
  for (auto i = 0; i < model()->columnCount(); ++i)
    cells.append(new QStandardItem(""));

  parentItem->appendRow(cells);
  auto newRowIndex = model()->index(parentItem->rowCount() - 1, 0, parent);
  m_subscriber->notifyRowInserted(rowLocationAt(newRowIndex));
  return newRowIndex;
}

QModelIndex JobTreeView::appendedSiblingRow(QModelIndex const &index,
                                            QList<QStandardItem *> cells) {
  return appendedChildRow(model()->parent(index), cells);
}

QModelIndex JobTreeView::appendedEmptySiblingRow(QModelIndex const &index) {
  return appendedEmptyChildRow(model()->parent(index));
}

QModelIndex JobTreeView::insertedChildRow(QModelIndex const &parent, int row,
                                          QList<QStandardItem *> cells) {
  auto *const parentItem = modelItemFromIndex(parent);
  parentItem->insertRow(row, cells);
  return model()->index(row, 0, parent);
}

QModelIndex JobTreeView::expanded(QModelIndex const &index) {
  auto expandAt = index;
  while (expandAt.isValid()) {
    setExpanded(expandAt, true);
    expandAt = model()->parent(expandAt);
  }
  return index;
}

QModelIndex JobTreeView::editAt(QModelIndex const &index) {
  clearSelection();
  setCurrentIndex(index);
  edit(index);
  return index;
}

QtTreeCursorNavigation JobTreeView::navigation() const {
  return QtTreeCursorNavigation(model());
}

QModelIndex JobTreeView::findOrMakeCellBelow(QModelIndex const &index) {
  if (navigation().isNotLastRowInThisNode(index)) {
    return moveCursor(QAbstractItemView::MoveDown, Qt::NoModifier);
  } else {
    return appendedEmptySiblingRow(index);
  }
}

void JobTreeView::keyPressEvent(QKeyEvent *event) {
  if (event->key() == Qt::Key_Return) {
    event->accept();
    if (event->modifiers() & Qt::ControlModifier)
      editAt(appendedEmptyChildRow(currentIndex()));
    else {
      auto below = findOrMakeCellBelow(currentIndex());
      editAt(expanded(below));
    }
  } else if (event->key() == Qt::Key_E) {
    static auto filter = false;
    static auto *originalModel = model();
    filter = !filter;
    if (filter) {
      auto *proxy = new QtFilterLeafNodes();
      proxy->setDynamicSortFilter(true);
      proxy->setSourceModel(originalModel);
      setModel(proxy);
    } else
      setModel(originalModel);
  } else {
    QTreeView::keyPressEvent(event);
  }
}

QModelIndex
JobTreeView::applyNavigationResult(QtTreeCursorNavigationResult const &result) {
  auto shouldMakeNewRowBelow = result.first;
  if (shouldMakeNewRowBelow)
    return expanded(appendedEmptySiblingRow(result.second));
  else
    return result.second;
}

QModelIndex JobTreeView::moveCursor(CursorAction cursorAction,
                                    Qt::KeyboardModifiers modifiers) {
  if (cursorAction == QAbstractItemView::MoveNext) {
    return applyNavigationResult(navigation().moveCursorNext(currentIndex()));
  } else if (cursorAction == QAbstractItemView::MovePrevious) {
    return navigation().moveCursorPrevious(currentIndex());
  } else {
    return QTreeView::moveCursor(cursorAction, modifiers);
  }
}

} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
