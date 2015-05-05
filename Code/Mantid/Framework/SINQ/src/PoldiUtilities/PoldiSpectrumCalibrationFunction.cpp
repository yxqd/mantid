#include "MantidSINQ/PoldiUtilities/PoldiSpectrumCalibrationFunction.h"
#include "MantidAPI/FunctionFactory.h"

namespace Mantid {
namespace Poldi {

using namespace API;

DECLARE_FUNCTION(PoldiSpectrumCalibrationFunction)

PoldiSpectrumCalibrationFunction::PoldiSpectrumCalibrationFunction()
    : PoldiSpectrumDomainFunction() {}

void PoldiSpectrumCalibrationFunction::init() {
  setDecoratedFunction("PoldiCalibrationProfile");
}

void PoldiSpectrumCalibrationFunction::functionModificationHook(
    const Poldi2DHelper_sptr &poldi2DHelper) const {
  m_profileFunction->setAttribute(
      "DeltaTheta", IFunction::Attribute(poldi2DHelper->deltaTwoTheta / 2.0));
}

} // namespace Poldi
} // namespace Mantid
