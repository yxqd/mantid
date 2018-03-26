#ifndef MANTIDQTMANTIDWIDGETS_QTTREEMODELADAPTER_H_
#define MANTIDQTMANTIDWIDGETS_QTTREEMODELADAPTER_H_
#include "QtTreeRow.h"
#include "TreeNode.h"
#include <QAbstractItemModel>
#include <QStringList>

namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

using Cell = QtTreeCell;
using RowNode = TreeNode<QtTreeRow>;

class QtTreeModelAdapter : public QAbstractItemModel {
  Q_OBJECT
public:
  QtTreeModelAdapter(QObject *parent = nullptr);
  QModelIndex index(int row, int column,
                    const QModelIndex &parent = QModelIndex()) const override;
  QModelIndex parent(const QModelIndex &index) const override;
  int rowCount(const QModelIndex &parent = QModelIndex()) const override;
  int columnCount(const QModelIndex &parent = QModelIndex()) const override;
  QVariant data(const QModelIndex &index,
                int role = Qt::DisplayRole) const override;
  QVariant headerData(int section, Qt::Orientation orientation,
                      int role) const override;
  bool setData(const QModelIndex &index, const QVariant &value,
               int role = Qt::EditRole) override;
  Qt::ItemFlags flags(const QModelIndex &index) const override;
  QModelIndex insertRow(const QModelIndex &parent, QtTreeRow row);

  virtual ~QtTreeModelAdapter() = default;

  Cell &cellAt(const QModelIndex &modelIndex);

private:
  RowNode &
  nodeFromModelIndex(const QModelIndex &modelIndex = QModelIndex()) const;

  RowNode *nthChildIfExists(RowNode &parent, std::size_t childIndex) const;

  bool isLeaf(QModelIndex const &modelIndex) const;

  mutable RowNode m_root;
};
} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
#endif // MANTIDQTMANTIDWIDGETS_QTTREEMODELADAPTER_H_
