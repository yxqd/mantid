#include "MantidQtWidgets/Common/Batch/QtTreeModelAdapter.h"
#include <algorithm>
#include <iostream>

namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

QtTreeModelAdapter::QtTreeModelAdapter(QObject *parent)
    : QAbstractItemModel(parent),
      m_root({Cell("Apples"), Cell("Bananas")},
             {RowNode({Cell("A"), Cell("B")},
                      {RowNode({Cell("H"), Cell("I")},
                               {RowNode({Cell("J"), Cell("K")})})}),
              RowNode({Cell("E"), Cell("F")})}) {
  m_root.updateParentPointers();
}

RowNode &
QtTreeModelAdapter::nodeFromModelIndex(const QModelIndex &modelIndex) const {
  if (!modelIndex.isValid()) // Invalid => parent is the root.
    return m_root;
  else
    return *static_cast<RowNode *>(modelIndex.internalPointer());
}

RowNode *QtTreeModelAdapter::nthChildIfExists(RowNode &parent,
                                              std::size_t childIndex) const {
  if (childIndex < parent.size())
    return &parent[childIndex];
  else
    return nullptr;
}

QModelIndex QtTreeModelAdapter::index(int row, int column,
                                      const QModelIndex &parentIndex) const {

  if (!hasIndex(row, column, parentIndex))
    return QModelIndex();

  auto &parentItem = nodeFromModelIndex(parentIndex);
  auto *childItem = nthChildIfExists(parentItem, row);
  if (childItem != nullptr)
    return createIndex(row, column, childItem);
  else
    return QModelIndex();
}

QModelIndex
QtTreeModelAdapter::parent(const QModelIndex &childModelIndex) const {
  auto &child = nodeFromModelIndex(childModelIndex);
  if (isRoot(child)) {
    return QModelIndex();
  } else {
    auto &parent = child.parent();
    if (isRoot(parent))
      return createIndex(0, 0, &parent);
    else
      return createIndex(indexOfNodeWithinParent(parent), 0, &parent);
  }
}

bool QtTreeModelAdapter::isLeaf(QModelIndex const &modelIndex) const {
  return modelIndex.column() == 0;
}

int QtTreeModelAdapter::rowCount(const QModelIndex &parent) const {
  return static_cast<int>(nodeFromModelIndex(parent).size());
}

int QtTreeModelAdapter::columnCount(const QModelIndex &) const { return 2; }

QVariant QtTreeModelAdapter::data(const QModelIndex &index, int role) const {
  if (index.isValid() && (role == Qt::DisplayRole || role == Qt::EditRole)) {
    auto &node = nodeFromModelIndex(index);
    auto &rowData = node.value();
    if (index.column() < rowData.size())
      return rowData[index.column()].text();
  }
  return QVariant();
}

QVariant QtTreeModelAdapter::headerData(int section,
                                        Qt::Orientation orientation,
                                        int role) const {
  if (orientation == Qt::Horizontal && role == Qt::DisplayRole)
    return m_root.value()[section].text();
  else
    return QVariant();
}

bool QtTreeModelAdapter::setData(const QModelIndex &index,
                                 const QVariant &value, int role) {
  if (role == Qt::EditRole) {
    auto &node = nodeFromModelIndex(index);
    auto &rowData = node.value();
    rowData[index.column()].setText(value.toString());
    emit dataChanged(index, index);
    return true;
  } else {
    return false;
  }
}

Qt::ItemFlags QtTreeModelAdapter::flags(const QModelIndex &index) const {
  if (!index.isValid())
    return 0;
  return Qt::ItemIsEditable | QAbstractItemModel::flags(index);
}

Cell& QtTreeModelAdapter::cellAt(const QModelIndex &modelIndex) {
  auto& row = nodeFromModelIndex(modelIndex);
  return row.value()[modelIndex.column()];
}

}; // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
