#include "MantidParallel/ThreadingBackend.h"

namespace Mantid {
namespace Parallel {
namespace detail {

ThreadingBackend::ThreadingBackend(const int size) : m_size(size) {}

int ThreadingBackend::size() const { return m_size; }

void ThreadingBackend::barrier() const {
  throw std::runtime_error(
      "ThreadingBackend::barrier() is not implemented yet.");
}

} // namespace detail
} // namespace Parallel
} // namespace Mantid
