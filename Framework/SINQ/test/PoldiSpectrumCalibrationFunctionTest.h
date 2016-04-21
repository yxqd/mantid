#ifndef MANTID_SINQ_POLDISPECTRUMCALIBRATIONFUNCTIONTEST_H_
#define MANTID_SINQ_POLDISPECTRUMCALIBRATIONFUNCTIONTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidSINQ/PoldiUtilities/PoldiSpectrumCalibrationFunction.h"

using Mantid::Poldi::PoldiSpectrumCalibrationFunction;
using namespace Mantid::API;

class PoldiSpectrumCalibrationFunctionTest : public CxxTest::TestSuite
{
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static PoldiSpectrumCalibrationFunctionTest *createSuite() { return new PoldiSpectrumCalibrationFunctionTest(); }
  static void destroySuite( PoldiSpectrumCalibrationFunctionTest *suite ) { delete suite; }


  void test_Something()
  {
    TSM_ASSERT( "You forgot to write a test!", 0);
  }


};


#endif /* MANTID_SINQ_POLDISPECTRUMCALIBRATIONFUNCTIONTEST_H_ */
