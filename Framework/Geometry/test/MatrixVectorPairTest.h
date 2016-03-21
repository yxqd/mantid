#ifndef MANTID_GEOMETRY_MATRIXVECTORPAIRTEST_H_
#define MANTID_GEOMETRY_MATRIXVECTORPAIRTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidGeometry/Crystal/MatrixVectorPair.h"
#include "MantidGeometry/Crystal/V3R.h"
#include "MantidKernel/Matrix.h"

using namespace Mantid::Geometry;
using namespace Mantid::Kernel;

typedef MatrixVectorPair<int, V3R> V3RIntPair;

class MatrixVectorPairTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static MatrixVectorPairTest *createSuite() {
    return new MatrixVectorPairTest();
  }
  static void destroySuite(MatrixVectorPairTest *suite) { delete suite; }

  void test_Construction() {
    TS_ASSERT_THROWS_NOTHING(V3RIntPair pair);
    TS_ASSERT_THROWS_NOTHING(
        V3RIntPair pair(Eigen::Matrix3i::Identity(), V3R()));
  }

  void test_getMatrix() {
    Eigen::Matrix3i m(Eigen::Matrix3i::Identity());
    m(0, 2) = 10;
    m(1, 1) = 5;
    m(2, 0) = 3;

    V3RIntPair pair(m, V3R());
    Eigen::Matrix3i inPair = pair.getMatrix();

    TS_ASSERT_EQUALS(inPair, m);
  }

  void test_getVector() {
    V3R v(2, 3, 4);

    V3RIntPair pair(Eigen::Matrix3i::Identity(), v);

    V3R inPair = pair.getVector();
    TS_ASSERT_EQUALS(inPair, v);
  }

  void test_multiplicationOperatorVector() {
    Eigen::Matrix3i m(Eigen::Matrix3i::Identity());
    m *= -1;
    V3R v(1, 1, 1);

    V3RIntPair pair(m, v);

    V3R toTransform(2, 3, 4);

    // m * t -> -t, t' + v -> t' + 1,1,1 -> -1, -2, -3
    V3R transformed = pair * toTransform;

    TS_ASSERT_EQUALS(transformed[0], -1);
    TS_ASSERT_EQUALS(transformed[1], -2);
    TS_ASSERT_EQUALS(transformed[2], -3);
  }

  void test_multiplicationOperatorMatrixVectorPair() {
    Eigen::Matrix3i m(Eigen::Matrix3i::Identity());
    m *= -1;

    V3RIntPair pairLHS(m, V3R(1, 1, 1));
    V3RIntPair pairRHS(m, V3R(2, 3, 4));

    V3RIntPair result = pairLHS * pairRHS;
    TS_ASSERT_EQUALS(result.getMatrix(), Eigen::Matrix3i::Identity());
    TS_ASSERT_EQUALS(result.getVector(), V3R(-1, -2, -3));
  }

  void test_inverse() {
    Eigen::Matrix3i m(Eigen::Matrix3i::Identity());
    m *= -1;
    m(1, 1) = 1;

    V3RIntPair pair(m, V3R(1, 2, 3));
    V3RIntPair inverse = pair.getInverse();

    V3R v(-2, 3, 1);
    TS_ASSERT_EQUALS(inverse * (pair * v), v);
  }

  void test_equalOperator() {
    Eigen::Matrix3i m(Eigen::Matrix3i::Identity());
    m *= -1;

    V3RIntPair pairLHS(m, V3R(1, 1, 1));
    V3RIntPair pairLHSCopy(pairLHS);

    V3RIntPair pairRHS(m, V3R(2, 3, 4));

    TS_ASSERT(pairLHS != pairRHS);
    TS_ASSERT_EQUALS(pairLHS, pairLHSCopy);
  }
};

#endif /* MANTID_GEOMETRY_MATRIXVECTORPAIRTEST_H_ */
