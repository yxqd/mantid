#ifndef MANTID_DATAOBJECTS_NULLMDBOXTEST_H_
#define MANTID_DATAOBJECTS_NULLMDBOXTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidDataObjects/NullMDBox.h"

using Mantid::DataObjects::NullMDBox;

class NullMDBoxTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static NullMDBoxTest *createSuite() { return new NullMDBoxTest(); }
  static void destroySuite( NullMDBoxTest *suite ) { delete suite; }


  void test_Something()
  {
    TS_FAIL( "You forgot to write a test!");
  }


};


#endif /* MANTID_DATAOBJECTS_NULLMDBOXTEST_H_ */