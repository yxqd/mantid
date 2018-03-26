#include "MantidQtWidgets/Common/Batch/QtTreeCell.h"
namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

QtTreeCell::QtTreeCell(QString text)
    : m_text(std::move(text)), m_borderThickness(2), m_borderColor(Qt::darkGray) {}

void QtTreeCell::setText(QString const &text) {
  m_text = text;
}

QString const &QtTreeCell::text() const {
  return m_text;
}

void QtTreeCell::setBorderThickness(int thickness) {
  m_borderThickness = thickness;
}

int QtTreeCell::borderThickness() const {
  return m_borderThickness;
}

void QtTreeCell::setBorderColor(QColor const &borderColor) {
  m_borderColor = borderColor;
}

QColor const &QtTreeCell::borderColor() const {
  return m_borderColor;
}

} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
