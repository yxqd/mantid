#ifndef MANTID_GEOMETRY_SPACEGROUPFACTORYTEST_H_
#define MANTID_GEOMETRY_SPACEGROUPFACTORYTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidGeometry/Crystal/SpaceGroupFactory.h"
#include "MantidGeometry/Crystal/SpaceGroupRegistration.h"

using namespace Mantid::Geometry;
using namespace Mantid::Geometry::SpaceGroupRegistration;

class SpaceGroupRegistrationTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static SpaceGroupRegistrationTest *createSuite() {
    return new SpaceGroupRegistrationTest();
  }
  static void destroySuite(SpaceGroupRegistrationTest *suite) { delete suite; }

  void test_SpaceGroupSubscriber_subscribes_to_factory() {
    TestableSpaceGroupFactory factory;

    GeneratedSpaceGroupSubscriber subscriber(factory);

    TS_ASSERT(!factory.isSubscribed("P 1"));
    TS_ASSERT(!factory.isSubscribed("P1"));

    subscriber(1, "P 1", "x,y,z");

    TS_ASSERT(factory.isSubscribed("P 1"));
    TS_ASSERT(factory.isSubscribed("P1"));

    TS_ASSERT_EQUALS(factory.subscribedSpaceGroupNumbers().size(), 1);
  }

  void test_SpaceGroupSubscriber_subscribes_aliases() {
    TestableSpaceGroupFactory factory;

    GeneratedSpaceGroupSubscriber subscriber(factory);

    TS_ASSERT(!factory.isSubscribed("test"));
    TS_ASSERT(!factory.isSubscribed("tests"));

    subscriber(1, "P 1", "x,y,z", "test,tests");

    TS_ASSERT(factory.isSubscribed("test"));
    TS_ASSERT(factory.isSubscribed("tests"));

    TS_ASSERT_EQUALS(factory.createSpaceGroup("test")->hmSymbol(), "P 1");
  }

  void test_SpaceGroupSubscriber_subscribes_aliases_without_spaces() {
    TestableSpaceGroupFactory factory;

    GeneratedSpaceGroupSubscriber subscriber(factory);

    TS_ASSERT(!factory.isSubscribed("test a"));
    TS_ASSERT(!factory.isSubscribed("testa"));

    subscriber(1, "P 1", "x,y,z", "test a");

    TS_ASSERT(factory.isSubscribed("test a"));
    TS_ASSERT(factory.isSubscribed("testa"));

    TS_ASSERT_EQUALS(factory.createSpaceGroup("test a")->hmSymbol(), "P 1");
  }

  void test_OrthorhombicSpaceGroupSubscriber_subscribes_to_factory() {
    TestableSpaceGroupFactory factory;

    OrthorhombicSpaceGroupSubscriber subscriber(factory);

    TS_ASSERT_EQUALS(factory.subscribedSpaceGroupNumbers().size(), 0);

    subscriber(17, "P 2 2 21", "-x,-y,z+1/2; -x,y,-z+1/2");

    TS_ASSERT_EQUALS(factory.subscribedSpaceGroupNumbers().size(), 1);
    TS_ASSERT_EQUALS(factory.subscribedSpaceGroupSymbols().size(), 3);

    TS_ASSERT(factory.isSubscribed("P 2 2 21"));
    TS_ASSERT(factory.isSubscribed("P2221"));
    TS_ASSERT(factory.isSubscribed("P 2 21 2"));
    TS_ASSERT(factory.isSubscribed("P2212"));
    TS_ASSERT(factory.isSubscribed("P 21 2 2"));
    TS_ASSERT(factory.isSubscribed("P2122"));
  }

private:
  class TestableSpaceGroupFactory : public SpaceGroupFactoryImpl {
    friend class SpaceGroupFactoryTest;

  public:
    TestableSpaceGroupFactory() : SpaceGroupFactoryImpl() {}
    ~TestableSpaceGroupFactory() override {}
  };
};

#endif /* MANTID_GEOMETRY_SPACEGROUPFACTORYTEST_H_ */
