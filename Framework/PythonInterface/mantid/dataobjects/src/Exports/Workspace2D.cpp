#include "MantidDataObjects/Workspace2D.h"
#include "MantidPythonInterface/kernel/Converters/NDArrayToVector.h"
#include "MantidPythonInterface/kernel/GetPointer.h"
#include "MantidPythonInterface/kernel/NdArray.h"
#include "MantidPythonInterface/kernel/Registry/RegisterWorkspacePtrToPython.h"

#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/raw_function.hpp>

using Mantid::API::MatrixWorkspace;
using Mantid::DataObjects::Workspace2D;
using Mantid::PythonInterface::Registry::RegisterWorkspacePtrToPython;
using Mantid::PythonInterface::Converters::NDArrayToVector;
namespace Np = Mantid::PythonInterface::NumPy;
namespace bp = boost::python;

GET_POINTER_SPECIALIZATION(Workspace2D)

namespace {
/// Alias for member function to modify the data
typedef Mantid::MantidVec &(MatrixWorkspace::*DataModifierFn)(
    const std::size_t);

/**
 * Set all of the 2D data in a single shot
 * @param self :: A reference to the calling object
 * @param accessor :: A member-function pointer to the data{X,Y,E} member that
 * will extract the writable values for a single spectrum
 * @param values :: A 2D numpy array. It's shape must match the shape of the
 * workspace.
 */
void set2DFromPyObject(MatrixWorkspace &self, DataModifierFn accessor,
                       const Np::NdArray &values) {
  if (values.get_nd() != 2) {
    throw std::invalid_argument("Expected a 2D numpy array, found nd=" +
                                std::to_string(values.get_nd()));
  }
  auto npshape = values.get_shape();
  if (static_cast<size_t>(npshape[0]) != self.getNumberHistograms()) {
    throw std::invalid_argument(
        "Expected first dimension to match number of histograms, found nd=" +
        std::to_string(npshape[0]));
  }
  for (size_t i = 0; i < self.getNumberHistograms(); ++i) {
    NDArrayToVector<double> converter(Np::NdArray(values[i]));
    converter.copyTo((self.*accessor)(i));
  }
}

/**
 * @param self  A reference to the calling object
 * @param x A 2D numpy array. It's shape must match the shape of the
 * workspace.
 * @param y A 2D numpy array. It's shape must match the shape of the
 * workspace.
 * @param e A 2D numpy array. It's shape must match the shape of the
 * workspace.
 * @param dx A 2D numpy array. It's shape must match the shape of the
 * workspace.
 */
void setDataOld(Workspace2D &self, const Np::NdArray &x = Np::NdArray(),
                const Np::NdArray &y = Np::NdArray(),
                const Np::NdArray &e = Np::NdArray(),
                const Np::NdArray &dx = Np::NdArray()) {
  size_t nonNullArgs(0);
  if (!x.is_none()) {
    set2DFromPyObject(self, &MatrixWorkspace::dataX, x);
    ++nonNullArgs;
  }
  if (!y.is_none()) {
    set2DFromPyObject(self, &MatrixWorkspace::dataY, y);
    ++nonNullArgs;
  }
  if (!e.is_none()) {
    set2DFromPyObject(self, &MatrixWorkspace::dataE, e);
    ++nonNullArgs;
  }
  if (!dx.is_none()) {
    set2DFromPyObject(self, &MatrixWorkspace::dataDx, dx);
    ++nonNullArgs;
  }
  if (nonNullArgs == 0) {
    throw std::invalid_argument(
        "At least 1 numpy array is required, none provided");
  }
}

/**
 * Emulates a python function defined as (*args, **kwargs). It does not
 * seem possible to attach a free C function as a member function in python
 * and still allow keyword arguments. The signature should be
 * setData(x=None, y=None, e=None, dx=None)
 * @param args Positional arguments
 * @param kwargs Keyword arguments
 */
bp::object setData(const bp::tuple &args, const bp::dict &kwargs) {
  // boost python ensures this is the case
  Workspace2D &self = bp::extract<Workspace2D &>(args[0]);
  auto nargs = bp::len(args) - 1;
  auto nkwargs = bp::len(kwargs);
  if (nargs == 0 && nkwargs == 0) {
    throw std::invalid_argument(
        "At least 1 numpy array is required, none provided");
  }
  static std::array<DataModifierFn, 4> methods{
      &MatrixWorkspace::dataX, &MatrixWorkspace::dataY, &MatrixWorkspace::dataE,
      &MatrixWorkspace::dataDx};
  if (nargs > 0) {
    for (size_t i = 0; i < methods.size(); ++i) {
      set2DFromPyObject(self, methods[i],
                        bp::extract<Np::NdArray>(args[i + 1]));
    }
  }

  //  switch (nargs) {
  //  case 1:
  //    set2DFromPyObject(self, &MatrixWorkspace::dataX,
  //                      bp::extract<Np::NdArray>(args[1]));
  //  case 2:
  //    set2DFromPyObject(self, &MatrixWorkspace::dataY,
  //                      bp::extract<Np::NdArray>(args[2]));

  //  case 4:
  //    // all positional arguments
  //    set2DFromPyObject(self, &MatrixWorkspace::dataX,
  //                      bp::extract<Np::NdArray>(args[1]));
  //    set2DFromPyObject(self, &MatrixWorkspace::dataY,
  //                      bp::extract<Np::NdArray>(args[2]));
  //    set2DFromPyObject(self, &MatrixWorkspace::dataE,
  //                      bp::extract<Np::NdArray>(args[3]));
  //    set2DFromPyObject(self, &MatrixWorkspace::dataDx,
  //                      bp::extract<Np::NdArray>(args[4]));
  //  }
  return bp::object(); // None
}
}

void export_Workspace2D() {
  bp::class_<Workspace2D, bp::bases<MatrixWorkspace>, boost::noncopyable>(
      "Workspace2D", bp::no_init)
      .def("setData", bp::raw_function(setData, 1),
           "Set X/Y/E/Dx values for a workspace");

  // register pointers
  RegisterWorkspacePtrToPython<Workspace2D>();
}
