#ifndef MANTID_MDALGORITHMS_RANKRESPONSIBILITYTEST_H_
#define MANTID_MDALGORITHMS_RANKRESPONSIBILITYTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidMDAlgorithms/RankResponsibility.h"

using Mantid::MDAlgorithms::RankResponsibility;

class RankResponsibilityTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static RankResponsibilityTest *createSuite() { return new RankResponsibilityTest(); }
  static void destroySuite( RankResponsibilityTest *suite ) { delete suite; }


  void test_Something()
  {
    TS_FAIL( "You forgot to write a test!");
  }


};


#endif /* MANTID_MDALGORITHMS_RANKRESPONSIBILITYTEST_H_ */