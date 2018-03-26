#ifndef MANTIDQTMANTIDWIDGETS_QTTREEROW_H_
#define MANTIDQTMANTIDWIDGETS_QTTREEROW_H_
#include "MantidQtWidgets/Common/DllOption.h"
#include "MantidQtWidgets/Common/Batch/QtTreeCell.h"
#include <QString>
#include <vector>
#include <initializer_list>

namespace MantidQt {
namespace MantidWidgets {
namespace Batch {
class EXPORT_OPT_MANTIDQT_COMMON QtTreeRow {
public:
  using iterator = typename std::vector<QtTreeCell>::iterator;
  using const_iterator = typename std::vector<QtTreeCell>::const_iterator;

  explicit QtTreeRow(std::vector<QtTreeCell> cells);
  QtTreeRow(std::initializer_list<QtTreeCell> const& cells);

  iterator begin();
  const_iterator cbegin() const;
  const_iterator begin() const;

  iterator end();
  const_iterator cend() const;
  const_iterator end() const;

  std::size_t size() const;

  QtTreeCell &operator[](std::size_t index);
  QtTreeCell const &operator[](std::size_t index) const;
private:
  std::vector<QtTreeCell> m_cells;
};
} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
#endif // MANTIDQTMANTIDWIDGETS_QTTREEROW_H_
