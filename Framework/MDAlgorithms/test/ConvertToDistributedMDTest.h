#ifndef MANTID_MDALGORITHMS_CONVERTTODISTRIBUTEDMDTEST_H_
#define MANTID_MDALGORITHMS_CONVERTTODISTRIBUTEDMDTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidMDAlgorithms/ConvertToDistributedMD.h"
#include "MantidTestHelpers/WorkspaceCreationHelper.h"


using Mantid::MDAlgorithms::ConvertToDistributedMD;

class ConvertToDistributedMDTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static ConvertToDistributedMDTest *createSuite() { return new ConvertToDistributedMDTest(); }
  static void destroySuite( ConvertToDistributedMDTest *suite ) { delete suite; }


  void test_Init()
  {
//    ConvertToDistributedMD alg;
//    TS_ASSERT_THROWS_NOTHING( alg.initialize() )
//    TS_ASSERT( alg.isInitialized() )
  }

  void test_exec() {
    // Arrange
    auto ws = generate_test_workspace(2);
    ConvertToDistributedMD alg;
    alg.initialize();
    alg.setProperty("InputWorkspace", ws);
    alg.setProperty("OutputWorkspace", "dummy"); // Not used for now, but have to set.

    // Act
    alg.execute();

    // Assert
    int a = 1;
  }

private:
  Mantid::DataObjects::EventWorkspace_sptr generate_test_workspace(int numberOfBanks) {
    // Create simple event workspace
    // This will create numberOfBanks pixels with 100 bins with 200 events. With pulse times
    // 0, 1, 2, ...
    auto ws = WorkspaceCreationHelper::createEventWorkspaceWithFullInstrument(numberOfBanks, 1, false /* clear events*/);
    return ws;
  }
};


#endif /* MANTID_MDALGORITHMS_CONVERTTODISTRIBUTEDMDTEST_H_ */