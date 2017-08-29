#include "MantidQtWidgets/MplCpp/NDArray1D.h"
#include "MantidQtWidgets/MplCpp/PythonErrors.h"
#include "MantidQtWidgets/Common/PythonThreading.h"

#include "MantidKernel/WarningSuppressions.h"

// See https://docs.scipy.org/doc/numpy/reference/c-api.array.html#miscellaneous
#define PY_ARRAY_UNIQUE_SYMBOL MPLCPP_ARRAY_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <algorithm>
#include <type_traits>
#include <vector>

namespace MantidQt {
namespace Widgets {
namespace MplCpp {

namespace {
// Simply struct to aid in global initialization of numpy
struct ImportArray {
  ImportArray() {
    PythonGIL gil;
    int result = _import_array();
    if (result < 0) {
      throw PythonError();
    }
  }
  ~ImportArray() {}
};
// static instance to call import_array()
// Do not remove this!!
// ImportArray _importer;

void initializeNumpy() {
  static ImportArray importer;
  (void)importer;
}
}

namespace detail {

// Numpy macro expands to code block containing a warning.
// clang-format off
GCC_DIAG_OFF(pedantic)
// clang-format on
template <typename Iterable> PyObject *copyToNDArray(const Iterable &data) {
  static_assert(std::is_same<typename Iterable::value_type, double>::value,
                "Element type must be double.");
  initializeNumpy();
  npy_intp length = static_cast<npy_intp>(data.size());
  auto ndarray = PyArray_SimpleNew(1, &length, NPY_DOUBLE);
  if (!ndarray)
    throw PythonError();
  auto emptyData =
      static_cast<double *>(PyArray_DATA((PyArrayObject *)ndarray));
  std::copy(std::begin(data), std::end(data), emptyData);
  return ndarray;
}
// clang-format off
GCC_DIAG_ON(pedantic)
// clang-format on

// Explicit template instantiations
template EXPORT_OPT_MANTIDQT_MPLCPP PyObject *
copyToNDArray<std::vector<double>>(const std::vector<double> &);
}

/**
 * Access the shape of the array
 * @return A single element array with the length of the array
 */
std::array<size_t, 1> NDArray1D::shape() const {
  auto npShape = PyArray_SHAPE((PyArrayObject *)(this->get()));
  return {{static_cast<size_t>(npShape[0])}};
}
}
}
}
