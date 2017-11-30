#ifndef MANTID_MDALGORITHMS_EVENTTOMDEVENTCONVERTERTEST_H_
#define MANTID_MDALGORITHMS_EVENTTOMDEVENTCONVERTERTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidMDAlgorithms/EventToMDEventConverter.h"
#include "MantidTestHelpers/WorkspaceCreationHelper.h"


using Mantid::MDAlgorithms::EventToMDEventConverter;

class EventToMDEventConverterTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static EventToMDEventConverterTest *createSuite() { return new EventToMDEventConverterTest(); }
  static void destroySuite( EventToMDEventConverterTest *suite ) { delete suite; }

  void test_that_can_convert_all_events_in_q_lab() {
    // Arrange
    auto ws = generate_test_workspace(2);
    EventToMDEventConverter converter;

    // Act
    auto mdEvents = converter.getEvents(*ws);

    // Assert
    // TODO: Test conversion correctness.
    TS_ASSERT(mdEvents.size() == 400);
  }

  void test_that_can_convert_half_of_the_events_in_q_lab() {
    // Arrange
    auto ws = generate_test_workspace(2);
    EventToMDEventConverter converter;

    // Act
    auto mdEvents = converter.getEvents(*ws, 0.5);

    // Assert
    // TODO: Test conversion correctness.
    TS_ASSERT(mdEvents.size() == 200);
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


#endif /* MANTID_MDALGORITHMS_EVENTTOMDEVENTCONVERTERTEST_H_ */