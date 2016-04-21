#include "MantidCurveFitting/Functions/PoldiCalibrationProfile.h"
#include "MantidAPI/FunctionFactory.h"
#include <iostream>

namespace Mantid {
namespace CurveFitting {
namespace Functions {

using namespace API;

DECLARE_FUNCTION(PoldiCalibrationProfile)

/**
 * Calculates function values
 *
 * This function calculates the function values like a Gaussian, but with some
 * additional modifications, as detailed in the class description.
 *
 * @param out :: The calculated function values.
 * @param xValues :: The x-values for which to calculate the function.
 * @param nData :: Number of data points.
 */
void PoldiCalibrationProfile::functionLocal(double *out, const double *xValues,
                                            const size_t nData) const {

  const double height = getParameter("Height");
  const double peakCentre = getParameter("PeakCentre");
  const double shift = getShiftFactor();
  const double realCentre =
      peakCentre * (1.0 + shift) + getAttribute("ChopperOffset").asDouble();

  const double weight = pow(1 / getParameter("Sigma"), 2);

  for (size_t i = 0; i < nData; i++) {
    double diff = xValues[i] - realCentre;
    out[i] = height * exp(-0.5 * diff * diff * weight);
  }
}

/**
 * Calculates derivatives
 *
 * This function calculates the Jacobian like a gaussian, but with some
 * additional modifications, as detailed in the class description.
 *
 * @param out :: The calculated Jacobian.
 * @param xValues :: The x-values for which to calculate the function.
 * @param nData :: Number of data points.
 */
void PoldiCalibrationProfile::functionDerivLocal(Jacobian *out,
                                                 const double *xValues,
                                                 const size_t nData) {
  const double height = getParameter("Height");
  const double peakCentre = getParameter("PeakCentre");
  const double shift = getShiftFactor();
  const double realCentre =
      peakCentre * (1.0 + shift) + getAttribute("ChopperOffset").asDouble();

  const double weight = pow(1 / getParameter("Sigma"), 2);
  const double factor = getConstantFactor();

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

/// Calulates the shift factor by multiplying the slope paramater by the return
/// value of getConstantFactor.
double PoldiCalibrationProfile::getShiftFactor() const {
  return getParameter("Slope") * getConstantFactor();
}

/// Returns the constant factor for the slope, depends on the DeltaTheta
/// attribute.
double PoldiCalibrationProfile::getConstantFactor() const {
  return getAttribute("DeltaTheta").asDouble() * 1.e-3;
}

/// Initialize Gaussian parameters and declare additional parameter.
void PoldiCalibrationProfile::init() {
  Gaussian::init();

  declareParameter("Slope", 0.0);
  declareAttribute("DeltaTheta", IFunction::Attribute(0.0));
  declareAttribute("ChopperOffset", IFunction::Attribute(0.0));
}

}

} // namespace CurveFitting
} // namespace Mantid
