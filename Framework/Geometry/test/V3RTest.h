#ifndef MANTID_GEOMETRY_V3RTEST_H_
#define MANTID_GEOMETRY_V3RTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidGeometry/Crystal/V3R.h"
#include "MantidKernel/Exception.h"

using namespace Mantid::Geometry;
using Mantid::Kernel::V3D;
using Mantid::Kernel::IntMatrix;

#include "Eigen/Core"
#include "MantidKernel/MersenneTwister.h"

class V3RTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static V3RTest *createSuite() { return new V3RTest(); }
  static void destroySuite(V3RTest *suite) { delete suite; }

  void testConstructors() {
    // default constructor
    V3R defConstr;
    TS_ASSERT_EQUALS(defConstr.x(), 0);
    TS_ASSERT_EQUALS(defConstr.y(), 0);
    TS_ASSERT_EQUALS(defConstr.z(), 0);

    // Constructor from rational numbers
    V3R rational(RationalNumber(1, 4), RationalNumber(1, 2),
                 RationalNumber(2, 3));
    V3D rationalV3D = rational;
    TS_ASSERT_EQUALS(rationalV3D.X(), 0.25);
    TS_ASSERT_EQUALS(rationalV3D.Y(), 0.5);
    TS_ASSERT_EQUALS(rationalV3D.Z(), 2.0 / 3.0);

    std::vector<int> good(3, 1);
    V3R rationalIntVec(good);
    TS_ASSERT_EQUALS(rationalIntVec.x(), 1);
    TS_ASSERT_EQUALS(rationalIntVec.y(), 1);
    TS_ASSERT_EQUALS(rationalIntVec.z(), 1);

    std::vector<int> bad(4, 1);
    TS_ASSERT_THROWS(V3R rationalIntVecBad(bad),
                     Mantid::Kernel::Exception::MisMatch<size_t>);

    // copy constructor
    V3R copied(rational);
    TS_ASSERT_EQUALS(copied.x(), rational.x());
    TS_ASSERT_EQUALS(copied.y(), rational.y());
    TS_ASSERT_EQUALS(copied.z(), rational.z());
  }

  void testXGetterSetter() {
    V3R vector;
    TS_ASSERT_EQUALS(vector.x(), 0);

    TS_ASSERT_THROWS_NOTHING(vector.setX(RationalNumber(1, 4)));
    TS_ASSERT_EQUALS(vector.x(), RationalNumber(1, 4));
  }

  void testYGetterSetter() {
    V3R vector;
    TS_ASSERT_EQUALS(vector.y(), 0);

    TS_ASSERT_THROWS_NOTHING(vector.setY(RationalNumber(1, 4)));
    TS_ASSERT_EQUALS(vector.y(), RationalNumber(1, 4));
  }

  void testZGetterSetter() {
    V3R vector;
    TS_ASSERT_EQUALS(vector.z(), 0);

    TS_ASSERT_THROWS_NOTHING(vector.setZ(RationalNumber(1, 4)));
    TS_ASSERT_EQUALS(vector.z(), RationalNumber(1, 4));
  }

  void testArrayAccess() {
    V3R vector(1, 2, 3);
    TS_ASSERT_EQUALS(vector[0], 1);
    TS_ASSERT_EQUALS(vector[1], 2);
    TS_ASSERT_EQUALS(vector[2], 3);
    TS_ASSERT_THROWS(vector[3], Mantid::Kernel::Exception::IndexError);

    TS_ASSERT_THROWS_NOTHING(vector[0] = RationalNumber(2, 3));
    TS_ASSERT_THROWS_NOTHING(vector[1] = RationalNumber(2, 3));
    TS_ASSERT_THROWS_NOTHING(vector[2] = RationalNumber(2, 3));
    TS_ASSERT_THROWS(vector[3] = RationalNumber(2, 3),
                     Mantid::Kernel::Exception::IndexError);
  }

  void testIntegerAddition() {
    V3R vector(RationalNumber(1, 4), RationalNumber(2, 3),
               RationalNumber(1, 2));
    V3R originalVector(vector);

    V3R vectorAdd = vector + 1;
    TS_ASSERT_EQUALS(vectorAdd.x(), RationalNumber(5, 4));
    TS_ASSERT_EQUALS(vectorAdd.y(), RationalNumber(5, 3));
    TS_ASSERT_EQUALS(vectorAdd.z(), RationalNumber(3, 2));

    vector += 1;
    TS_ASSERT_EQUALS(vector, vectorAdd);

    vector += -1;
    TS_ASSERT_EQUALS(vector, originalVector);
  }

  void testIntegerSubtraction() {
    V3R vector(RationalNumber(1, 4), RationalNumber(2, 3),
               RationalNumber(1, 2));
    V3R originalVector(vector);

    V3R vectorAdd = vector - 1;
    TS_ASSERT_EQUALS(vectorAdd.x(), RationalNumber(-3, 4));
    TS_ASSERT_EQUALS(vectorAdd.y(), RationalNumber(-1, 3));
    TS_ASSERT_EQUALS(vectorAdd.z(), RationalNumber(-1, 2));

    vector -= 1;
    TS_ASSERT_EQUALS(vector, vectorAdd);

    vector -= -1;
    TS_ASSERT_EQUALS(vector, originalVector);
  }

  void testIntegerMultiplication() {
    V3R vector(RationalNumber(1, 4), RationalNumber(2, 3),
               RationalNumber(1, 2));
    V3R originalVector(vector);

    V3R vectorAdd = vector * 2;
    TS_ASSERT_EQUALS(vectorAdd.x(), RationalNumber(1, 2));
    TS_ASSERT_EQUALS(vectorAdd.y(), RationalNumber(4, 3));
    TS_ASSERT_EQUALS(vectorAdd.z(), RationalNumber(1, 1));

    vector *= 2;
    TS_ASSERT_EQUALS(vector, vectorAdd);

    vector /= 2;
    TS_ASSERT_EQUALS(vector, originalVector);

    V3R nullVector = vector * 0;
    TS_ASSERT_EQUALS(nullVector.x(), 0);
    TS_ASSERT_EQUALS(nullVector.y(), 0);
    TS_ASSERT_EQUALS(nullVector.z(), 0);
  }

  void testIntegerDivision() {
    V3R vector(RationalNumber(1, 4), RationalNumber(2, 3),
               RationalNumber(1, 2));
    V3R originalVector(vector);

    V3R vectorAdd = vector / 2;
    TS_ASSERT_EQUALS(vectorAdd.x(), RationalNumber(1, 8));
    TS_ASSERT_EQUALS(vectorAdd.y(), RationalNumber(1, 3));
    TS_ASSERT_EQUALS(vectorAdd.z(), RationalNumber(1, 4));

    vector /= 2;
    TS_ASSERT_EQUALS(vector, vectorAdd);

    vector *= 2;
    TS_ASSERT_EQUALS(vector, originalVector);

    TS_ASSERT_THROWS(vector / 0, boost::bad_rational);
  }

  void testRationalAddition() {
    V3R vector(RationalNumber(1, 4), RationalNumber(2, 3),
               RationalNumber(1, 2));
    V3R originalVector(vector);

    V3R vectorAdd = vector + RationalNumber(1, 2);
    TS_ASSERT_EQUALS(vectorAdd.x(), RationalNumber(3, 4));
    TS_ASSERT_EQUALS(vectorAdd.y(), RationalNumber(7, 6));
    TS_ASSERT_EQUALS(vectorAdd.z(), RationalNumber(1, 1));

    vector += RationalNumber(1, 2);
    TS_ASSERT_EQUALS(vector, vectorAdd);

    vector += RationalNumber(-1, 2);
    TS_ASSERT_EQUALS(vector, originalVector);
  }

  void testRationalSubtraction() {
    V3R vector(RationalNumber(1, 4), RationalNumber(2, 3),
               RationalNumber(1, 2));
    V3R originalVector(vector);

    V3R vectorAdd = vector - RationalNumber(1, 2);
    TS_ASSERT_EQUALS(vectorAdd.x(), RationalNumber(-1, 4));
    TS_ASSERT_EQUALS(vectorAdd.y(), RationalNumber(1, 6));
    TS_ASSERT_EQUALS(vectorAdd.z(), RationalNumber(0));

    vector -= RationalNumber(1, 2);
    TS_ASSERT_EQUALS(vector, vectorAdd);

    vector -= RationalNumber(-1, 2);
    TS_ASSERT_EQUALS(vector, originalVector);
  }

  void testRationalMultiplication() {
    V3R vector(RationalNumber(1, 4), RationalNumber(2, 3),
               RationalNumber(1, 2));
    V3R originalVector(vector);

    V3R vectorAdd = vector * RationalNumber(1, 2);
    TS_ASSERT_EQUALS(vectorAdd.x(), RationalNumber(1, 8));
    TS_ASSERT_EQUALS(vectorAdd.y(), RationalNumber(1, 3));
    TS_ASSERT_EQUALS(vectorAdd.z(), RationalNumber(1, 4));

    vector *= RationalNumber(1, 2);
    TS_ASSERT_EQUALS(vector, vectorAdd);

    vector /= RationalNumber(1, 2);
    TS_ASSERT_EQUALS(vector, originalVector);
  }

  void testRationalDivision() {
    V3R vector(RationalNumber(1, 4), RationalNumber(2, 3),
               RationalNumber(1, 2));
    V3R originalVector(vector);

    V3R vectorAdd = vector / RationalNumber(1, 2);
    TS_ASSERT_EQUALS(vectorAdd.x(), RationalNumber(1, 2));
    TS_ASSERT_EQUALS(vectorAdd.y(), RationalNumber(4, 3));
    TS_ASSERT_EQUALS(vectorAdd.z(), RationalNumber(1));

    vector /= RationalNumber(1, 2);
    TS_ASSERT_EQUALS(vector, vectorAdd);

    vector *= RationalNumber(1, 2);
    TS_ASSERT_EQUALS(vector, originalVector);
  }

  void testVectorAddition() {
    V3R vector(RationalNumber(1, 4), RationalNumber(2, 3),
               RationalNumber(1, 2));
    V3R otherVector(RationalNumber(-3, 7), RationalNumber(1, 3),
                    RationalNumber(7, 9));
    V3R originalVector(vector);

    V3R vectorAdd = vector + otherVector;
    TS_ASSERT_EQUALS(vectorAdd.x(), RationalNumber(-5, 28));
    TS_ASSERT_EQUALS(vectorAdd.y(), RationalNumber(1));
    TS_ASSERT_EQUALS(vectorAdd.z(), RationalNumber(23, 18));

    vector += otherVector;
    TS_ASSERT_EQUALS(vector, vectorAdd);

    vector += -otherVector;
    TS_ASSERT_EQUALS(vector, originalVector);

    V3R nullVector = vector + (-vector);
    TS_ASSERT_EQUALS(nullVector.x(), 0);
    TS_ASSERT_EQUALS(nullVector.y(), 0);
    TS_ASSERT_EQUALS(nullVector.z(), 0);
  }

  void testVectorSubtraction() {
    V3R vector(RationalNumber(1, 4), RationalNumber(2, 3),
               RationalNumber(1, 2));
    V3R otherVector(RationalNumber(-3, 7), RationalNumber(1, 3),
                    RationalNumber(7, 9));
    V3R originalVector(vector);

    V3R vectorAdd = vector - otherVector;
    TS_ASSERT_EQUALS(vectorAdd.x(), RationalNumber(19, 28));
    TS_ASSERT_EQUALS(vectorAdd.y(), RationalNumber(1, 3));
    TS_ASSERT_EQUALS(vectorAdd.z(), RationalNumber(-5, 18));

    vector -= otherVector;
    TS_ASSERT_EQUALS(vector, vectorAdd);

    vector -= -otherVector;
    TS_ASSERT_EQUALS(vector, originalVector);

    V3R nullVector = vector - vector;
    TS_ASSERT_EQUALS(nullVector.x(), 0);
    TS_ASSERT_EQUALS(nullVector.y(), 0);
    TS_ASSERT_EQUALS(nullVector.z(), 0);
  }

  void testV3DAddition() {
    V3R vector(RationalNumber(1, 4), RationalNumber(2, 3),
               RationalNumber(1, 2));
    V3D factor(0.5, 0.5, 0.5);

    V3D newVector = factor + vector;

    TS_ASSERT_EQUALS(newVector.X(), 0.75);

    // It's not exactly equal because of floating point precision
    TS_ASSERT_DIFFERS(newVector.Y(), 7.0 / 6.0);
    TS_ASSERT_EQUALS(newVector.Y(), 0.5 + 2.0 / 3.0);
    TS_ASSERT_DELTA(newVector.Y(), 7.0 / 6.0, 1e-15);

    TS_ASSERT_EQUALS(newVector.Z(), 1.0);

    // check operation with different operand ordering
    V3D equalVector = vector + factor;
    TS_ASSERT_EQUALS(equalVector, newVector);
  }

  void testV3DSubtraction() {
    V3R vector(RationalNumber(1, 4), RationalNumber(2, 3),
               RationalNumber(1, 2));
    V3D factor(0.5, 0.5, 0.5);

    V3D newVector = factor - vector;

    TS_ASSERT_EQUALS(newVector.X(), 0.25);

    // It's not exactly equal because of floating point precision
    TS_ASSERT_DIFFERS(newVector.Y(), -1.0 / 6.0);
    TS_ASSERT_EQUALS(newVector.Y(), 0.5 - 2.0 / 3.0);
    TS_ASSERT_DELTA(newVector.Y(), -1.0 / 6.0, 1e-16);

    TS_ASSERT_EQUALS(newVector.Z(), 0.0);
  }

  void testEqualityOperator() {
    V3R one(RationalNumber(1, 4), RationalNumber(2, 3), RationalNumber(1, 2));
    V3R two(RationalNumber(1, 4), RationalNumber(2, 3), RationalNumber(1, 2));
    TS_ASSERT_EQUALS(one, two);

    V3R three(RationalNumber(2, 8), RationalNumber(6, 9),
              RationalNumber(14, 28));
    TS_ASSERT_EQUALS(one, three);

    V3R four(RationalNumber(1, 5), RationalNumber(2, 3), RationalNumber(1, 2));
    TS_ASSERT_DIFFERS(one, four);

    V3R five(RationalNumber(1, 4), RationalNumber(2, 4), RationalNumber(1, 2));
    TS_ASSERT_DIFFERS(one, five);

    V3R six(RationalNumber(1, 4), RationalNumber(2, 3), RationalNumber(1, 3));
    TS_ASSERT_DIFFERS(one, six);
  }

  void testComparison() {
    V3R one(RationalNumber(1, 4), RationalNumber(2, 3), RationalNumber(1, 2));

    V3R two(RationalNumber(1, 5), RationalNumber(2, 3), RationalNumber(1, 2));
    TS_ASSERT_LESS_THAN(two, one);

    V3R three(RationalNumber(1, 3), RationalNumber(2, 3), RationalNumber(1, 2));
    TS_ASSERT_LESS_THAN(one, three);

    V3R four(RationalNumber(1, 4), RationalNumber(2, 4), RationalNumber(1, 2));
    TS_ASSERT_LESS_THAN(four, one);

    V3R five(RationalNumber(1, 4), RationalNumber(2, 2), RationalNumber(1, 2));
    TS_ASSERT_LESS_THAN(one, five);

    V3R six(RationalNumber(1, 4), RationalNumber(2, 3), RationalNumber(1, 3));
    TS_ASSERT_LESS_THAN(six, one);

    V3R seven(RationalNumber(1, 4), RationalNumber(2, 3), RationalNumber(2, 2));
    TS_ASSERT_LESS_THAN(one, seven);
  }

  void testIntegerComparison() {
    V3R zeros;
    TS_ASSERT_EQUALS(zeros, 0);
    zeros.setX(RationalNumber(1, 2));
    TS_ASSERT_DIFFERS(zeros, 0);
  }

  void testMatrixMultiplication() {
    V3R vector(RationalNumber(1, 4), RationalNumber(2, 3),
               RationalNumber(1, 2));
    // unit matrix - resulting vector must be equal
    IntMatrix unity(3, 3, true);
    V3R transformedUnity = unity * vector;
    TS_ASSERT_EQUALS(transformedUnity, vector);

    // inversion
    IntMatrix inversion = unity * -1;
    V3R transformedInversion = inversion * vector;
    TS_ASSERT_EQUALS(transformedInversion, -vector);

    // general
    IntMatrix operation(3, 3);
    operation[0][0] = 0;
    operation[0][1] = 1;
    operation[0][2] = 1;

    operation[1][0] = 1;
    operation[1][1] = -1;
    operation[1][2] = 1;

    operation[2][0] = -1;
    operation[2][1] = -1;
    operation[2][2] = 0;

    V3R transformedGeneral = operation * vector;
    TS_ASSERT_EQUALS(transformedGeneral.x(), RationalNumber(7, 6)); // y + z
    TS_ASSERT_EQUALS(transformedGeneral.y(),
                     RationalNumber(1, 12)); // x - y + z
    TS_ASSERT_EQUALS(transformedGeneral.z(), RationalNumber(-11, 12)); // -x - y

    // wrong sizes
    IntMatrix wrongOne(3, 4);
    TS_ASSERT_THROWS(wrongOne * vector,
                     Mantid::Kernel::Exception::MisMatch<size_t>);

    IntMatrix wrongTwo(4, 3);
    TS_ASSERT_THROWS(wrongTwo * vector, Mantid::Kernel::Exception::IndexError);

    // Smaller works
    IntMatrix wrongThree(2, 3);
    wrongThree[0][0] = 1;
    wrongThree[0][1] = 0;
    wrongThree[0][2] = 0;

    wrongThree[1][0] = 0;
    wrongThree[1][1] = 1;
    wrongThree[1][2] = 0;

    TS_ASSERT_THROWS_NOTHING(wrongThree * vector);
    V3R transformedSmaller = wrongThree * vector;

    TS_ASSERT_EQUALS(transformedSmaller.x(), vector.x());
    TS_ASSERT_EQUALS(transformedSmaller.y(), vector.y());
    TS_ASSERT_EQUALS(transformedSmaller.z(), 0);
  }

  void testVectorOperator() {
    V3R test = V3R(1, 2, 3) / 4;

    std::vector<double> approximations(test);

    TS_ASSERT_EQUALS(approximations.size(), 3);
    TS_ASSERT_EQUALS(approximations[0], 0.25);
    TS_ASSERT_EQUALS(approximations[1], 0.5);
    TS_ASSERT_EQUALS(approximations[2], 0.75);
  }
};

namespace Eigen {
template <> struct NumTraits<RationalNumber> : NumTraits<int> {
  typedef int Real;
  typedef double NonInteger;
  typedef RationalNumber Nested;

  enum {
    IsComplex = 0,
    IsInteger = 0,
    ReadCost = 2,
    AddCost = 200,
    MulCost = 200,
    IsSigned = 1,
    RequireInitialization = 0
  };
};

namespace internal {
template <typename T> struct scalar_product_traits<boost::rational<T>, T> {
  enum { Defined = 1 };
  typedef boost::rational<T> ReturnType;
};

template <typename T> struct scalar_product_traits<T, boost::rational<T>> {
  enum { Defined = 1 };
  typedef boost::rational<T> ReturnType;
};
}
}

typedef Eigen::Matrix<RationalNumber, 3, 3> Matrix3r;
typedef Eigen::Matrix<RationalNumber, 3, 1> Vector3r;

using namespace Mantid::Geometry;
using namespace Mantid::Kernel;

class V3RTestPerformance : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static V3RTestPerformance *createSuite() { return new V3RTestPerformance(); }
  static void destroySuite(V3RTestPerformance *suite) { delete suite; }

  V3RTestPerformance() {
    setupDataMantid();
    setupDataEigen();
  }

  void test_mult_mantid() {
    std::vector<V3R> results(num_matrices);

    for (size_t i = 0; i < num_iterations; ++i) {
      std::transform(
          matrices_mantid.cbegin(), matrices_mantid.cend(),
          vectors_mantid.cbegin(), results.begin(),
          [](const IntMatrix &lhs, const V3R &rhs) { return lhs * rhs; });

      TS_ASSERT_EQUALS(results.size(), matrices_mantid.size());
    }
  }

  void test_mult_eigen() {
    std::vector<Vector3r> results(num_matrices);

    for (size_t i = 0; i < num_iterations; ++i) {
      std::transform(matrices_eigen.cbegin(), matrices_eigen.cend(),
                     vectors_eigen.cbegin(), results.begin(),
                     [](const Eigen::Matrix3i &lhs, const Vector3r &rhs) {
                       return lhs * rhs;
                     });

      TS_ASSERT_EQUALS(results.size(), matrices_eigen.size());
    }
  }

private:
  void setupDataMantid() {
    using Mantid::Kernel::MersenneTwister;
    MersenneTwister rng(4000, -1.0, 1.0);

    std::fill_n(
        std::back_inserter(vectors_mantid), num_matrices,
        V3R(RationalNumber(2, 3), RationalNumber(1, 3), RationalNumber(1, 4)));

    std::generate_n(std::back_inserter(matrices_mantid), num_matrices,
                    [&rng]() {
                      std::vector<int> init;
                      std::generate_n(std::back_inserter(init), 9, [&]() {
                        return static_cast<int>(rng.nextValue() * 100.0);
                      });
                      return init;
                    });
  }

  void setupDataEigen() {
    std::transform(matrices_mantid.cbegin(), matrices_mantid.cend(),
                   std::back_inserter(matrices_eigen),
                   [](const IntMatrix &matrix) {
                     Eigen::Matrix3i m;
                     for (const auto &val : matrix.getVector()) {
                       m << val;
                     }

                     return m;
                   });
    std::fill_n(std::back_inserter(vectors_eigen), matrices_eigen.size(),
                Vector3r(RationalNumber(2, 3), RationalNumber(1, 3),
                         RationalNumber(1, 4)));
  }

  std::vector<IntMatrix> matrices_mantid;
  std::vector<V3R> vectors_mantid;

  std::vector<Eigen::Matrix3i> matrices_eigen;
  std::vector<Vector3r> vectors_eigen;

  size_t num_matrices{100};
  size_t num_iterations{10000};
};

#endif /* MANTID_GEOMETRY_V3RTEST_H_ */
