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

void PoldiSpectrumCalibrationFunction::functionModificationPreHook(
    const Poldi2DHelper_sptr &poldi2DHelper) const {
  m_profileFunction->setAttribute(
      "DeltaTheta", IFunction::Attribute(poldi2DHelper->deltaTwoTheta / 2.0));
}

void PoldiSpectrumCalibrationFunction::functionModificationPostHook(
    const Poldi2DHelper_sptr &poldi2DHelper) const {
  UNUSED_ARG(poldi2DHelper);

  m_profileFunction->setAttribute("DeltaTheta", IFunction::Attribute(0.0));
}

void
PoldiSpectrumCalibrationFunction::setPeakCenter(double newCenter,
                                                double chopperOffset) const {
  m_profileFunction->setCentre(newCenter);
  m_profileFunction->setAttribute("ChopperOffset",
                                  IFunction::Attribute(chopperOffset));
}

} // namespace Poldi
} // namespace Mantid
