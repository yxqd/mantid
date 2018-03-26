#include "MantidQtWidgets/Common/Batch/QtTreeRow.h"
namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

QtTreeRow::QtTreeRow(std::vector<QtTreeCell> cells)
    : m_cells(std::move(cells)) {}

QtTreeRow::QtTreeRow(std::initializer_list<QtTreeCell> const &cells)
    : m_cells(cells) {}

auto QtTreeRow::begin() -> iterator { return m_cells.begin(); }

auto QtTreeRow::cbegin() const -> const_iterator { return m_cells.cbegin(); }

auto QtTreeRow::begin() const -> const_iterator { return m_cells.begin(); }

auto QtTreeRow::end() -> iterator { return m_cells.end(); }

auto QtTreeRow::cend() const -> const_iterator { return m_cells.cend(); }

auto QtTreeRow::end() const -> const_iterator { return m_cells.end(); }

std::size_t QtTreeRow::size() const { return m_cells.size(); }

QtTreeCell &QtTreeRow::operator[](std::size_t index) { return m_cells[index]; }

QtTreeCell const &QtTreeRow::operator[](std::size_t index) const {
  return m_cells[index];
}
} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
