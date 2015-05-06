#include "MantidCurveFitting/PoldiCalibrationProfile.h"
#include "MantidAPI/FunctionFactory.h"
#include <iostream>

namespace Mantid {
namespace CurveFitting {

using namespace API;

DECLARE_FUNCTION(PoldiCalibrationProfile)

void PoldiCalibrationProfile::functionLocal(double *out, const double *xValues,
                                            const size_t nData) const {

  const double height = getParameter("Height");
  const double peakCentre = getParameter("PeakCentre");
  const double shift = getAbsoluteShift();
  const double realCentre =
      peakCentre * (1.0 + shift) + getAttribute("ChopperOffset").asDouble();

  const double weight = pow(1 / getParameter("Sigma"), 2);

  for (size_t i = 0; i < nData; i++) {
    double diff = xValues[i] - realCentre;
    out[i] = height * exp(-0.5 * diff * diff * weight);
  }
}

void PoldiCalibrationProfile::functionDerivLocal(Jacobian *out,
                                                 const double *xValues,
                                                 const size_t nData) {
  const double height = getParameter("Height");
  const double peakCentre = getParameter("PeakCentre");
  const double shift = getAbsoluteShift();
  const double realCentre =
      peakCentre * (1.0 + shift) + getAttribute("ChopperOffset").asDouble();

  const double weight = pow(1 / getParameter("Sigma"), 2);
  const double factor = getAttribute("DeltaTheta").asDouble() * 1.e-3;

  for (size_t i = 0; i < nData; i++) {
    double diff = xValues[i] - realCentre;
    double e = exp(-0.5 * diff * diff * weight);
    double b = e * height * diff * weight;
    out->set(i, 0, e);
    out->set(i, 1, b * (shift + 1.0));
    out->set(i, 2, -0.5 * diff * diff * height *
                       e); // derivative with respect to weight not sigma
    out->set(i, 3, b * factor * peakCentre);
  }
}

double PoldiCalibrationProfile::getAbsoluteShift() const {
  return getParameter("Slope") * getAttribute("DeltaTheta").asDouble() * 1.e-3;
}

/// Initialize Gaussian parameters and declare additional parameter.
void PoldiCalibrationProfile::init() {
  Gaussian::init();

  declareParameter("Slope", 0.0);
  declareAttribute("DeltaTheta", IFunction::Attribute(0.0));
  declareAttribute("ChopperOffset", IFunction::Attribute(0.0));

  // fix(3);
}

} // namespace CurveFitting
} // namespace Mantid
