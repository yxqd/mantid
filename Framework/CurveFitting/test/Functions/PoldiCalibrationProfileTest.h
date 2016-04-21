#ifndef MANTID_CURVEFITTING_POLDICALIBRATIONFUNCTIONTEST_H_
#define MANTID_CURVEFITTING_POLDICALIBRATIONFUNCTIONTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidCurveFitting/Jacobian.h"
#include "MantidCurveFitting/Functions/PoldiCalibrationProfile.h"

#include <boost/lexical_cast.hpp>

using namespace Mantid::CurveFitting;
using namespace Mantid::CurveFitting::Functions;
using namespace Mantid::API;

class PoldiCalibrationProfileTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static PoldiCalibrationProfileTest *createSuite() {
    return new PoldiCalibrationProfileTest();
  }
  static void destroySuite(PoldiCalibrationProfileTest *suite) { delete suite; }

  void testParameters() {
    PoldiCalibrationProfile fn;
    fn.initialize();

    Gaussian gaussian;
    gaussian.initialize();

    TS_ASSERT_EQUALS(fn.nParams(), gaussian.nParams() + 1);
  }

  void testFunctionValue() {
    PoldiCalibrationProfile fn;
    fn.initialize();
    fn.setParameter("PeakCentre", 1.25);
    fn.setParameter("Sigma", 0.01);
    fn.setParameter("Height", 5.0);
    fn.setParameter("Slope", 0.1);

    fn.setAttribute("DeltaTheta", IFunction::Attribute(0.1));
    fn.setAttribute("ChopperOffset", IFunction::Attribute(0.001));

    /* The parameter combination corresponds to a Gaussian with the following
     * peak centre:
     *    c = 1.25 * (1.0 + 0.1 * 0.1 * 1e-3) + 0.001 = 1.2510125
     */
    Gaussian gaussian;
    gaussian.initialize();
    gaussian.setParameter("PeakCentre", 1.2510125);
    gaussian.setParameter("Sigma", 0.01);
    gaussian.setParameter("Height", 5.0);

    FunctionDomain1DVector domain(1.1, 1.4, 100);
    FunctionValues valuesPoldi(domain);
    FunctionValues valuesGauss(domain);

    fn.function(domain, valuesPoldi);
    gaussian.function(domain, valuesGauss);

    for (size_t i = 0; i < valuesPoldi.size(); ++i) {
      TSM_ASSERT_DELTA("Problem in point " +
                           boost::lexical_cast<std::string>(i) + " at x=" +
                           boost::lexical_cast<std::string>(domain[i]) + ".",
                       valuesPoldi[i], valuesGauss[i], 1e-13);
    }
  }
};

#endif /* MANTID_CURVEFITTING_POLDICALIBRATIONFUNCTIONTEST_H_ */
