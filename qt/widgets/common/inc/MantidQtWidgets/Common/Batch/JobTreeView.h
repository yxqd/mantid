#ifndef MANTIDQTMANTIDWIDGETS_JOBTREEVIEW_H_
#define MANTIDQTMANTIDWIDGETS_JOBTREEVIEW_H_
#include "MantidQtWidgets/Common/DllOption.h"
#include <QTreeView>
#include <QStyledItemDelegate>
#include <QStandardItemModel>
#include <QPainter>
#include <iostream>
namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

class CellDelegate : public QStyledItemDelegate {
public:
  explicit CellDelegate(QObject *parent, QTreeView const &view)
      : QStyledItemDelegate(parent), m_view(view){};

  void paint(QPainter *painter, const QStyleOptionViewItem &option,
             const QModelIndex &index) const override {
    QStyledItemDelegate::paint(painter, option, index);
    painter->save();
    auto pen =
        (m_view.currentIndex() == index) ? QPen(Qt::black) : QPen(Qt::darkGray);
    pen.setWidth(1);
    painter->setPen(pen);
    painter->drawRect(option.rect.adjusted(1, 1, -1, -1));
    painter->restore();
  }

private:
  QTreeView const &m_view;
};

using RowPath = std::vector<int>;

class RowLocation {
public:
  RowLocation(RowPath path) : m_path(std::move(path)) {}
  RowPath const &path() const { return m_path; }
  int rowRelativeToParent() const { return m_path.back(); }
  bool isRoot() const { return m_path.empty(); };

private:
  RowPath m_path;
};

class IJobTreeViewSubscriber {
  void notifyCellChanged(RowLocation itemIndex, int column,
                         std::string newValue);
  void notifyRequestRowInserted(RowLocation itemIndex);
  void notifyRequestRowRemoved(RowLocation itemIndex);
  void notifySelectedRowsChanged(std::vector<RowLocation> const &selection);
};

inline void assertOrThrow(bool condition, std::string const &message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

class EXPORT_OPT_MANTIDQT_COMMON JobTreeView : public QTreeView {
  Q_OBJECT
public:
  QList<QStandardItem *>
  rowFromStringVector(std::vector<std::string> const &rowText) const;
  JobTreeView(QWidget *parent = nullptr);
  JobTreeView(QStringList const &columnHeadings, QWidget *parent = nullptr);

  void subscribe(IJobTreeViewSubscriber &notifyee);

  void insertChildRowOf(RowLocation const &parent, int beforeRow,
                        std::vector<std::string> const &rowText);
  void appendChildRowOf(RowLocation const &parent);
  void appendChildRowOf(RowLocation const &parentLocation,
                        std::vector<std::string> const &rowText);

  void removeRowAt(RowLocation const &location);

  std::vector<std::string> rowTextAt(RowLocation const &location) const;
  void setRowTextAt(RowLocation const &location,
                    std::vector<std::string> const &rowText);

  std::string textAt(RowLocation location, int column);
  void setTextAt(RowLocation location, int column, std::string const &cellText);

  QModelIndex moveCursor(CursorAction cursorAction,
                         Qt::KeyboardModifiers modifiers) override;

protected:
  void keyPressEvent(QKeyEvent *event) override;
  QModelIndex addExpandAndMoveToSibling(QModelIndex const &index);
  QModelIndex addExpandAndMoveToChild(QModelIndex const &index);
  QModelIndex addSibling(QModelIndex const &index);
  QModelIndex addChild(QModelIndex const &index);
  void expand(QModelIndex const &index);
  void expandAndMoveTo(QModelIndex const &index);

private:
  QModelIndex parentIndexOf(RowLocation const &location) const;
  QStandardItem *parentItemOf(RowLocation const &location) const;

  QModelIndex modelIndexAt(RowLocation const &location) const;
  QStandardItem *modelItemAt(RowLocation const &location) const;

  QStandardItem *modelItemFromIndex(QModelIndex const &location) const;

  QStandardItem *emptyRow() const;
  IJobTreeViewSubscriber *m_notifyee;
  QStandardItemModel m_model;
};
} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
#endif // MANTIDQTMANTIDWIDGETS_JOBTREEVIEW_H_
