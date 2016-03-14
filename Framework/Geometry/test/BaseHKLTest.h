#ifndef MANTID_GEOMETRY_BASEHKLTEST_H_
#define MANTID_GEOMETRY_BASEHKLTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidGeometry/Crystal/BaseHKL.h"

using namespace Mantid::Geometry;

class BaseHKLTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static BaseHKLTest *createSuite() { return new BaseHKLTest(); }
  static void destroySuite(BaseHKLTest *suite) { delete suite; }

  void test_Something() {
    ProHKL hkl(3.0, 3.0, 3.0);
    TS_ASSERT_EQUALS(hkl.h(), 3.0);
  }

  void test_addition() {
    FractionalHKL hklF(3.5, 2.5, 1.5);
    ProHKL one(1, 2, 3);
    FractionalHKL two = static_cast<FractionalHKL>(one) + hklF;

    TS_ASSERT_EQUALS(two, FractionalHKL(4.5, 4.5, 4.5));
  }

  void test_conversion() {
    ProHKL one(-1.3, 1.5, 0.5);
    IntegerHKL two(one);

    TS_ASSERT_EQUALS(two, IntegerHKL(-1, 2, 1));

    ProHKL diff = one - static_cast<ProHKL>(two);

    TS_ASSERT_EQUALS(diff, ProHKL(-0.3, -0.5, -0.5));
  }
};

#endif /* MANTID_GEOMETRY_BASEHKLTEST_H_ */
