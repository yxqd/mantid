#ifndef MANTIDQTMANTIDWIDGETS_BATCHWIDGET_H_
#define MANTIDQTMANTIDWIDGETS_BATCHWIDGET_H_
#include "MantidQtWidgets/Common/Batch/QtTreeModelAdapter.h"
#include "MantidQtWidgets/Common/DllOption.h"
#include "ui_BatchWidget.h"
namespace MantidQt {
namespace MantidWidgets {
namespace Batch {
class EXPORT_OPT_MANTIDQT_COMMON BatchWidget : public QWidget {
  Q_OBJECT
public:
  BatchWidget(QWidget *parent = nullptr) : QWidget(parent) {
    m_ui.setupUi(this);
  }

private:
  Ui::BatchWidget m_ui;
};
} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
#endif // MANTIDQTMANTIDWIDGETS_BATCHWIDGET_H_
