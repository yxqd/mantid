#ifndef MANTIDQTMANTIDWIDGETS_JOBTREEVIEW_H_
#define MANTIDQTMANTIDWIDGETS_JOBTREEVIEW_H_
#include "MantidQtWidgets/Common/Batch/QtTreeModelAdapter.h"
#include "MantidQtWidgets/Common/DllOption.h"
#include "QtTreeRow.h"
#include <QTreeView>
#include <QStyledItemDelegate>
#include <QPainter>
namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

class CellDelegate : public QStyledItemDelegate {
public:
  explicit CellDelegate(QObject *parent, QtTreeModelAdapter &modelAdapter)
      : QStyledItemDelegate(parent), m_modelAdapter(modelAdapter) {};

  void paint(QPainter *painter, const QStyleOptionViewItem &option,
             const QModelIndex &index) const override {
    QStyledItemDelegate::paint(painter, option, index);
    auto& cell = m_modelAdapter.cellAt(index);
    painter->save();
    auto pen = QPen(cell.borderColor());
    pen.setWidth(cell.borderThickness());
    painter->setPen(pen);
    painter->drawRect(option.rect);
    painter->restore();
  }

private:
  QtTreeModelAdapter &m_modelAdapter;
};

class EXPORT_OPT_MANTIDQT_COMMON JobTreeView : public QTreeView {
  Q_OBJECT
public:
  JobTreeView(QWidget *parent = nullptr);
  QModelIndex moveCursor(CursorAction cursorAction,
                         Qt::KeyboardModifiers modifiers) override;
protected:

private:
  QtTreeModelAdapter m_modelAdapter;
};
} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
#endif // MANTIDQTMANTIDWIDGETS_JOBTREEVIEW_H_
