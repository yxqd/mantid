#ifndef MANTIDQTMANTIDWIDGETS_BATCHWIDGET_H_
#define MANTIDQTMANTIDWIDGETS_BATCHWIDGET_H_
#include "ui_BatchWidget.h"
#include "MantidQtWidgets/Common/DllOption.h"
namespace MantidQt {
namespace MantidWidgets {
namespace Batch {
class QtBatchModelAdaptor : public QAbstractItemModel {
  Q_OBJECT
public:
  QtBatchModelAdaptor(QObject *parent)
      : QAbstractItemModel(
            parent) /*, items({"Root", "1stChild", "2ndChild"}),
        parent({-1, 0, 0})*/ {}

  QModelIndex index(int row, int column,
                    const QModelIndex &parent = QModelIndex()) const override {
    // if(hasIndex(row, column, parent)) {
    // int* parent = nullptr;
    // if (parent.isValid()) { // Invalid => parent is the root.
    // parent = &items[0];
    //} else {
    // parent = ((1 + row) < 3) ? &items[0] : nullptr; // Hacked traversal for
    // nth child if exists.
    //}
    //} else
    return QModelIndex();
  }

  QModelIndex parent(const QModelIndex &index) const override {
    // if (index.isValid()) {
    // auto *itemAtIndex = static_cast<int *>(index.internalPointer());
    // return itemAtIndex;
    //// return invalid index for root.
    //} else {
    return QModelIndex();
    //}
  }

  int rowCount(const QModelIndex &parent = QModelIndex()) const override {
    return 0;
  }

  int columnCount(const QModelIndex &parent = QModelIndex()) const override {
    return 1;
  }

  QVariant data(const QModelIndex &index,
                int role = Qt::DisplayRole) const override {
    return "";
  }

  bool setData(const QModelIndex &index, const QVariant &value,
               int role = Qt::EditRole) override {
    return false;
  }

  Qt::ItemFlags flags(const QModelIndex &index) const override {
    return Qt::ItemIsEditable | Qt::ItemIsSelectable;
  }

  virtual ~QtBatchModelAdaptor() = default;
};

class EXPORT_OPT_MANTIDQT_COMMON BatchWidget : QObject {
  Q_OBJECT
  BatchWidget() : m_adaptor(nullptr) { m_ui.viewTable->setModel(&m_adaptor); }

private:
  QtBatchModelAdaptor m_adaptor;
  Ui::BatchWidget m_ui;
};
}
}
}
#endif // MANTIDQTMANTIDWIDGETS_BATCHWIDGET_H_
