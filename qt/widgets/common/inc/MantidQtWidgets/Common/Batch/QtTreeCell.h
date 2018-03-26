#ifndef MANTIDQTMANTIDWIDGETS_QTTREECELL_H_
#define MANTIDQTMANTIDWIDGETS_QTTREECELL_H_
#include "MantidQtWidgets/Common/DllOption.h"
#include <QString>
#include <QColor>

namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

class EXPORT_OPT_MANTIDQT_COMMON QtTreeCell {
public:
  explicit QtTreeCell(QString text);

  void setText(QString const& text);
  QString const& text() const;

  void setBorderThickness(int thickness);
  int borderThickness() const;

  void setBorderColor(QColor const& thickness);
  QColor const& borderColor() const;

private:
  QString m_text;
  QColor m_borderColor;
  int m_borderThickness;
};
} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
#endif // MANTIDQTMANTIDWIDGETS_QTTREECELL_H_
