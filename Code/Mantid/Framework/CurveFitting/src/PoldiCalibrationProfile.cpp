#include "MantidCurveFitting/PoldiCalibrationProfile.h"
#include "MantidAPI/FunctionFactory.h"

namespace Mantid {
namespace CurveFitting {

using namespace API;

DECLARE_FUNCTION(PoldiCalibrationProfile)

void PoldiCalibrationProfile::functionLocal(double *out, const double *xValues,
                                            const size_t nData) const {

  const double height = getParameter("Height");
  const double peakCentre = getParameter("PeakCentre");
  const double shift = getAbsoluteShift();
  const double realCentre = peakCentre + peakCentre * shift;

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
  const double realCentre = peakCentre + peakCentre * shift;

  const double weight = pow(1 / getParameter("Sigma"), 2);
  const double factor = getAttribute("DeltaTheta").asDouble() / 1000.0;

  for (size_t i = 0; i < nData; i++) {
    double diff = xValues[i] - realCentre;
    double e = exp(-0.5 * diff * diff * weight);
    double b = e * height * diff * weight;
    out->set(i, 0, e);
    out->set(i, 1, b * (1.0 + shift));
    out->set(i, 2, -0.5 * diff * diff * height *
                       e); // derivative with respect to weight not sigma
    out->set(i, 3, b * factor * peakCentre);
  }
}

double PoldiCalibrationProfile::centre() const {
  return getParameter("PeakCentre") +
         getParameter("PeakCentre") * getAbsoluteShift();
}

double PoldiCalibrationProfile::getAbsoluteShift() const {
  return getParameter("Slope") * 1.e-3 * getAttribute("DeltaTheta").asDouble();
}

/// Initialize Gaussian parameters and declare additional parameter.
void PoldiCalibrationProfile::init() {
  Gaussian::init();

  declareParameter("Slope");
  declareAttribute("DeltaTheta", IFunction::Attribute(0.0));
}

} // namespace CurveFitting
} // namespace Mantid
