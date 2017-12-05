#ifndef MANTID_MDALGORITHMS_BOXSTRUCTURESERIALIZERTEST_H_
#define MANTID_MDALGORITHMS_BOXSTRUCTURESERIALIZERTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidMDAlgorithms/BoxStructureSerializer.h"
#include "MantidTestHelpers/MDEventsTestHelper.h"


using Mantid::MDAlgorithms::BoxStructureSerializer;
using namespace Mantid::DataObjects;

class BoxStructureSerializerTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static BoxStructureSerializerTest *createSuite() { return new BoxStructureSerializerTest(); }
  static void destroySuite( BoxStructureSerializerTest *suite ) { delete suite; }


  void test_that_can_serialize_box_structure() {
    // Arrange
    auto b = MDEventsTestHelper::makeMDBox3();
    auto g = new MDGridBox<MDLeanEvent<3>, 3>(b);

    // Act
    BoxStructureSerializer serializer;
    auto serialInformation = serializer.createSerializedBoxStructure(g, *g->getBoxController());

    // Assert
    TS_ASSERT_EQUALS(serialInformation.boxType.size(), 101ul)
    TS_ASSERT_EQUALS(serialInformation.boxType[0], 2)
    TS_ASSERT_EQUALS(serialInformation.boxType[1], 1)

    // Clean up
    const auto bcc = b->getBoxController();
    delete bcc;
    delete g;
    delete b;
  }


};


#endif /* MANTID_MDALGORITHMS_BOXSTRUCTURESERIALIZERTEST_H_ */