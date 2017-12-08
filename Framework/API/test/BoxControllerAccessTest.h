#ifndef MANTID_API_BOXCONTROLLERACCESSTEST_H_
#define MANTID_API_BOXCONTROLLERACCESSTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidAPI/BoxControllerAccess.h"

using Mantid::API::BoxControllerAccess;

class BoxControllerAccessTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static BoxControllerAccessTest *createSuite() { return new BoxControllerAccessTest(); }
  static void destroySuite( BoxControllerAccessTest *suite ) { delete suite; }


  void test_Something()
  {
    TS_FAIL( "You forgot to write a test!");
  }


};


#endif /* MANTID_API_BOXCONTROLLERACCESSTEST_H_ */