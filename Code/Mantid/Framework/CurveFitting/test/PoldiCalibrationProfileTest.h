#ifndef MANTID_CURVEFITTING_POLDICALIBRATIONFUNCTIONTEST_H_
#define MANTID_CURVEFITTING_POLDICALIBRATIONFUNCTIONTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidCurveFitting/Jacobian.h"
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
};

#endif /* MANTID_CURVEFITTING_POLDICALIBRATIONFUNCTIONTEST_H_ */
