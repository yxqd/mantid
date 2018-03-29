#ifndef MANTIDQTMANTIDWIDGETS_JOBTREEVIEW_H_
#define MANTIDQTMANTIDWIDGETS_JOBTREEVIEW_H_
#include "MantidQtWidgets/Common/DllOption.h"
#include "MantidQtWidgets/Common/Batch/QtTreeCursorNavigation.h"
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
    pen.setWidth((m_view.currentIndex() == index) ? 2 : 1);
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
public:
  virtual void notifyCellChanged(RowLocation itemIndex, int column,
                                 std::string newValue) = 0;
  virtual void notifyRowInserted(RowLocation itemIndex) = 0;
  virtual void notifyRowRemoved(RowLocation itemIndex) = 0;
  virtual void
  notifySelectedRowsChanged(std::vector<RowLocation> const &selection) = 0;
};

inline void assertOrThrow(bool condition, std::string const &message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

class EXPORT_OPT_MANTIDQT_COMMON JobTreeView : public QTreeView {
  Q_OBJECT
public:
  JobTreeView(QWidget *parent = nullptr);
  JobTreeView(QStringList const &columnHeadings, QWidget *parent = nullptr);

  void subscribe(IJobTreeViewSubscriber &subscriber);

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
  void make(QModelIndex const &){};

protected:
  void keyPressEvent(QKeyEvent *event) override;

private:
  QModelIndex appendedEmptySiblingRow(QModelIndex const &index);
  QModelIndex appendedSiblingRow(QModelIndex const &index,
                                 QList<QStandardItem *> cells);
  QModelIndex appendedEmptyChildRow(QModelIndex const &parent);
  QModelIndex appendedChildRow(QModelIndex const &parent,
                               QList<QStandardItem *> cells);
  QModelIndex insertedChildRow(QModelIndex const &parent, int column,
                               QList<QStandardItem *> cells);

  QModelIndex expanded(QModelIndex const &index);
  QModelIndex editAt(QModelIndex const &index);

  QtTreeCursorNavigation navigation() const;

  QModelIndex applyNavigationResult(QtTreeCursorNavigationResult const &result);
  QModelIndex findOrMakeCellBelow(QModelIndex const &index);

  QList<QStandardItem *>
  rowFromRowText(std::vector<std::string> const &rowText) const;
  std::vector<std::string> rowTextFromRow(QModelIndex firstCellIndex) const;

  QModelIndex parentIndexOf(RowLocation const &location, int column = 0) const;
  QStandardItem *parentItemOf(RowLocation const &location,
                              int column = 0) const;

  QModelIndex modelIndexAt(RowLocation const &location, int column = 0) const;
  RowLocation rowLocationAt(QModelIndex const &index) const;
  QStandardItem *modelItemAt(RowLocation const &location, int column = 0) const;

  QStandardItem *modelItemFromIndex(QModelIndex const &location) const;

  IJobTreeViewSubscriber *m_notifyee;
  QStandardItemModel m_model;
};
} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
#endif // MANTIDQTMANTIDWIDGETS_JOBTREEVIEW_H_
