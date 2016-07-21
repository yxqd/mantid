#include "MantidDataObjects/WorkspaceSingleValue.h"
#include "MantidAPI/WorkspaceProperty.h"
#include "MantidAPI/WorkspaceFactory.h"

namespace Mantid {
namespace DataObjects {

using std::size_t;

DECLARE_WORKSPACE(WorkspaceSingleValue)

/// Constructor
WorkspaceSingleValue::WorkspaceSingleValue(double value, double error)
    : API::MatrixWorkspace() {
  // Set the "histogram" to the single value
  data.dataX().resize(1, 0.0);
  data.setFrequencies(1, value);
  data.setFrequencyStandardDeviations(1, error);
  data.setPointStandardDeviations(1, 0.0);
}

/** Does nothing in this case
*  @param NVectors :: This value can only be equal to one, otherwise exception
* is thrown
*  @param XLength :: The number of X data points/bin boundaries
*  @param YLength :: The number of data/error points
*/
void WorkspaceSingleValue::init(const HistogramData::Histogram::YMode ymode,
                                const std::size_t &NVectors,
                                const std::size_t &XLength,
                                const std::size_t &YLength) {
  // There is no bin width for a single point so only Frequencies makes sense.
  if(ymode != HistogramData::Histogram::YMode::Frequencies)
    throw std::invalid_argument(
        "WorkspaceSingleValue can only be initialized with YMode::Frequencies");

  (void)NVectors;
  (void)XLength;
  (void)YLength; // Avoid compiler warning
}

/// Return the underlying Histogram1D at the given workspace index.
Histogram1D &WorkspaceSingleValue::getSpectrum(const size_t /*index*/) {
  return data;
}

/// Return the underlying Histogram1D at the given workspace index.
const Histogram1D &
WorkspaceSingleValue::getSpectrum(const size_t /*index*/) const {
  return data;
}

/// Rebin the workspace. Not implemented for this workspace.
void WorkspaceSingleValue::generateHistogram(const std::size_t index,
                                             const MantidVec &X, MantidVec &Y,
                                             MantidVec &E,
                                             bool skipError) const {
  UNUSED_ARG(index);
  UNUSED_ARG(X);
  UNUSED_ARG(Y);
  UNUSED_ARG(E);
  UNUSED_ARG(skipError);
  throw std::runtime_error(
      "generateHistogram() not implemented for WorkspaceSingleValue.");
}

/// Our parent MatrixWorkspace has hardcoded 2, but we need 0.
size_t WorkspaceSingleValue::getNumDims() const { return 0; }

} // namespace DataObjects
} // namespace Mantid

///\cond TEMPLATE

template class DLLExport
    Mantid::API::WorkspaceProperty<Mantid::DataObjects::WorkspaceSingleValue>;

namespace Mantid {
namespace Kernel {
template <>
DLLExport Mantid::DataObjects::WorkspaceSingleValue_sptr
IPropertyManager::getValue<Mantid::DataObjects::WorkspaceSingleValue_sptr>(
    const std::string &name) const {
  PropertyWithValue<Mantid::DataObjects::WorkspaceSingleValue_sptr> *prop =
      dynamic_cast<
          PropertyWithValue<Mantid::DataObjects::WorkspaceSingleValue_sptr> *>(
          getPointerToProperty(name));
  if (prop) {
    return *prop;
  } else {
    std::string message =
        "Attempt to assign property " + name +
        " to incorrect type. Expected shared_ptr<WorkspaceSingleValue>.";
    throw std::runtime_error(message);
  }
}

template <>
DLLExport Mantid::DataObjects::WorkspaceSingleValue_const_sptr
IPropertyManager::getValue<
    Mantid::DataObjects::WorkspaceSingleValue_const_sptr>(
    const std::string &name) const {
  PropertyWithValue<Mantid::DataObjects::WorkspaceSingleValue_sptr> *prop =
      dynamic_cast<
          PropertyWithValue<Mantid::DataObjects::WorkspaceSingleValue_sptr> *>(
          getPointerToProperty(name));
  if (prop) {
    return prop->operator()();
  } else {
    std::string message =
        "Attempt to assign property " + name +
        " to incorrect type. Expected const shared_ptr<WorkspaceSingleValue>.";
    throw std::runtime_error(message);
  }
}

} // namespace Kernel
} // namespace Mantid

///\endcond TEMPLATE
