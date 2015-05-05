#ifndef MANTID_CURVEFITTING_POLDICALIBRATIONFUNCTIONTEST_H_
#define MANTID_CURVEFITTING_POLDICALIBRATIONFUNCTIONTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidCurveFitting/PoldiCalibrationProfile.h"

using namespace Mantid::CurveFitting;
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

  void testDerivatives() {
      PoldiCalibrationProfile fn;
      fn.initialize();

      fn.setParameter("Height", 2.0);
      fn.setParameter("PeakCentre", 1.1);
      fn.setParameter("Sigma", 0.01);
      fn.setParameter("Slope", 0.001);
      fn.setAttribute("DeltaTheta", IFunction::Attribute(0.01));


  }
};

#endif /* MANTID_CURVEFITTING_POLDICALIBRATIONFUNCTIONTEST_H_ */
