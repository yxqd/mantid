#include "MantidDataObjects/EventWorkspaceMRU.h"
#include "MantidKernel/System.h"

namespace Mantid {
namespace DataObjects {

EventWorkspaceMRU::~EventWorkspaceMRU() {
  // Make sure you free up the memory in the MRUs
  {
    Poco::ScopedWriteRWLock _lock(m_changeMruListsMutexY);
    m_bufferedDataY.clear();
  }
  Poco::ScopedWriteRWLock _lock2(m_changeMruListsMutexE);
  m_bufferedDataE.clear();
}

//---------------------------------------------------------------------------
/// Clear all the data in the MRU buffers
void EventWorkspaceMRU::clear() {
  {
    // Make sure you free up the memory in the MRUs
    Poco::ScopedWriteRWLock _lock(m_changeMruListsMutexY);
    m_bufferedDataY.clear();
  }
  Poco::ScopedWriteRWLock _lock(m_changeMruListsMutexE);
  m_bufferedDataE.clear();
}

//---------------------------------------------------------------------------
/** Find a Y histogram in the MRU
 *
 * @param thread_num :: number of the thread in which this is run
 * @param index :: index of the data to return
 * @return pointer to the TypeWithMarker that has the data; NULL if not found.
 */
Kernel::cow_ptr<HistogramData::HistogramY>
EventWorkspaceMRU::findY(const EventList *index) {
  Poco::ScopedReadRWLock _lock(m_changeMruListsMutexY);
  auto result =
      m_bufferedDataY.find(reinterpret_cast<const std::uintptr_t>(index));
  if (result)
    return result->m_data;
  return YType(nullptr);
}

/** Find a Y histogram in the MRU
 *
 * @param thread_num :: number of the thread in which this is run
 * @param index :: index of the data to return
 * @return pointer to the TypeWithMarker that has the data; NULL if not found.
 */
Kernel::cow_ptr<HistogramData::HistogramE>
EventWorkspaceMRU::findE(const EventList *index) {
  Poco::ScopedReadRWLock _lock(m_changeMruListsMutexE);
  auto result =
      m_bufferedDataE.find(reinterpret_cast<const std::uintptr_t>(index));
  if (result)
    return result->m_data;
  return EType(nullptr);
}

/** Insert a new histogram into the MRU
 *
 * @param thread_num :: thread being accessed
 * @param data :: the new data
 * @param index :: index of the data to insert
 */
void EventWorkspaceMRU::insertY(YType data, const EventList *index) {
  Poco::ScopedReadRWLock _lock(m_changeMruListsMutexY);
  auto yWithMarker =
      new TypeWithMarker<YType>(reinterpret_cast<const std::uintptr_t>(index));
  yWithMarker->m_data = std::move(data);
  auto oldData = m_bufferedDataY.insert(yWithMarker);
  // And clear up the memory of the old one, if it is dropping out.
  delete oldData;
}

/** Insert a new histogram into the MRU
 *
 * @param thread_num :: thread being accessed
 * @param data :: the new data
 * @param index :: index of the data to insert
 */
void EventWorkspaceMRU::insertE(EType data, const EventList *index) {
  Poco::ScopedReadRWLock _lock(m_changeMruListsMutexE);
  auto eWithMarker =
      new TypeWithMarker<EType>(reinterpret_cast<const std::uintptr_t>(index));
  eWithMarker->m_data = std::move(data);
  auto oldData = m_bufferedDataE.insert(eWithMarker);
  // And clear up the memory of the old one, if it is dropping out.
  delete oldData;
}

/** Delete any entries in the MRU at the given index
 *
 * @param index :: index to delete.
 */
void EventWorkspaceMRU::deleteIndex(const EventList *index) {
  {
    Poco::ScopedReadRWLock _lock1(m_changeMruListsMutexE);
    m_bufferedDataE.deleteIndex(reinterpret_cast<const std::uintptr_t>(index));
  }
  Poco::ScopedReadRWLock _lock2(m_changeMruListsMutexY);
  m_bufferedDataY.deleteIndex(reinterpret_cast<const std::uintptr_t>(index));
}

size_t EventWorkspaceMRU::MRUSize() const {
  return this->m_bufferedDataY.size();
}

} // namespace Mantid
} // namespace DataObjects
