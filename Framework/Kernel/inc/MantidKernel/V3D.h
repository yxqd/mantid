#ifndef MANTID_KERNEL_V3D_H_
#define MANTID_KERNEL_V3D_H_

#include <cmath>
#include <cfloat>
#include <complex>
#include <vector>
#include "MantidKernel/DllConfig.h"
#include "MantidKernel/Matrix.h"
#include "MantidKernel/Tolerance.h"
#include "MantidKernel/Exception.h"
#include <nexus/NeXusFile.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Mantid {
namespace Kernel {
/** @class V3D V3D.h Kernel\V3D.h

Class for 3D vectors.

@author Laurent C Chapon, ISIS, RAL
@date 09/10/2007

Copyright &copy; 2007-8 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
National Laboratory & European Spallation Source

This file is part of Mantid.

Mantid is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Mantid is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

File change history is stored at: <https://github.com/mantidproject/mantid>.
Code Documentation is available at: <http://doxygen.mantidproject.org>
*/
class MANTID_KERNEL_DLL V3D {
public:
  V3D() : m_vector(0, 0, 0) {}
  V3D(const double x, const double y, const double z) : m_vector(x, y, z) {}
  explicit V3D(const Eigen::Vector3d &vector) : m_vector(vector) {}

  /// Convenience method for sorting list of V3D objects based on magnitude
  static bool CompareMagnitude(const Kernel::V3D &v1, const Kernel::V3D &v2);

  // explicit conversion into vector
  operator std::vector<double>() const {
    return std::vector<double>{m_vector.data(),
                               m_vector.data() + m_vector.size()};
  }

  // Arithemetic operators overloaded
  V3D operator+(const V3D &v) const { return V3D(m_vector + v.m_vector); }

  V3D &operator+=(const V3D &v) {
    m_vector += v.m_vector;
    return *this;
  }

  V3D operator-(const V3D &v) const { return V3D(m_vector - v.m_vector); }

  V3D &operator-=(const V3D &v) {
    m_vector -= v.m_vector;
    return *this;
  }

  V3D operator*(const V3D &v) const {
    return V3D(m_vector.cwiseProduct(v.m_vector));
  }

  V3D operator/(const V3D &v) const {
    return V3D(m_vector.cwiseQuotient(v.m_vector));
  }

  V3D &operator*=(const V3D &v) {
    m_vector = m_vector.cwiseProduct(v.m_vector);
    return *this;
  }

  V3D &operator/=(const V3D &v) {
    if (!v.m_vector.isZero()) {
      m_vector = m_vector.cwiseQuotient(v.m_vector);
    }
    return *this;
  }

  // Scale
  V3D operator*(const double D) const { return V3D(m_vector * D); }

  V3D &operator*=(const double D) {
    m_vector *= D;
    return *this;
  }

  V3D operator/(const double D) const {
    V3D out(*this);
    out /= D;
    return out;
  }

  V3D &operator/=(const double D) {
    if (D != 0.0) {
      m_vector /= D;
    }
    return *this;
  }

  // Simple Comparison
  bool operator==(const V3D &v) const {
    return (m_vector - v.m_vector).isZero(Kernel::Tolerance);
  }

  bool operator!=(const V3D &other) const { return !(this->operator==(other)); }

  bool operator<(const V3D &) const;
  bool operator>(const V3D &rhs) const;

  // Access
  // Setting x, y and z values
  void operator()(const double x, const double y, const double z) {
    m_vector = Eigen::Vector3d(x, y, z);
  }

  void spherical(const double &R, const double &theta, const double &phi);
  void spherical_rad(const double &R, const double &polar,
                     const double &azimuth);
  void azimuth_polar_SNS(const double &R, const double &azimuth,
                         const double &polar);
  void setX(const double x) { m_vector(0) = x; }
  void setY(const double y) { m_vector(1) = y; }
  void setZ(const double z) { m_vector(2) = z; }

  const double &X() const { return m_vector(0); } ///< Get x
  const double &Y() const { return m_vector(1); } ///< Get y
  const double &Z() const { return m_vector(2); } ///< Get z

  const double &operator[](const size_t Index) const {
    if (Index > 2) {
      throw Kernel::Exception::IndexError(Index, 2,
                                          "V3D::operator[] range error");
    }
    return m_vector(Index);
  }

  double &operator[](const size_t Index) {
    if (Index > 2) {
      throw Kernel::Exception::IndexError(Index, 2,
                                          "V3D::operator[] range error");
    }
    return m_vector(Index);
  }

  void getSpherical(double &R, double &theta, double &phi) const;

  //      void rotate(const V3D&,const V3D&,const double);
  void rotate(const Matrix<double> &A)
  /**
  Rotate a point by a matrix
  @param A :: Rotation matrix (needs to be >3x3)
*/
  {
    double xold(m_vector(0)), yold(m_vector(1)), zold(m_vector(2));
    m_vector(0) = A[0][0] * xold + A[0][1] * yold + A[0][2] * zold;
    m_vector(1) = A[1][0] * xold + A[1][1] * yold + A[1][2] * zold;
    m_vector(2) = A[2][0] * xold + A[2][1] * yold + A[2][2] * zold;
  }

  void round() {
    m_vector = m_vector.unaryExpr(
        [](const double &element) { return std::round(element); });
  }

  /// Make a normalized vector (return norm value)
  double normalize() {
    const double ND(norm());
    this->operator/=(ND);
    return ND;
  } // Vec3D::makeUnit
  double norm() const { return m_vector.norm(); }
  double norm2() const { return m_vector.squaredNorm(); }
  /// transform vector into form, used to describe directions in
  /// crystallogaphical coodinate system
  double toMillerIndexes(double eps = 1.e-3);
  /// Scalar product
  double scalar_prod(const V3D &v) const { return m_vector.dot(v.m_vector); }
  /// Cross product
  V3D cross_prod(const V3D &v) const {
    V3D out(*this);
    out.m_vector = m_vector.cross(v.m_vector);
    return out;
  }
  /// Distance (R) between two points defined as vectors
  double distance(const V3D &v) const { return (m_vector - v.m_vector).norm(); }
  /// Zenith (theta) angle between this and another vector
  double zenith(const V3D &) const;
  /// Angle between this and another vector
  double angle(const V3D &) const;
  /// Direction angles
  V3D directionAngles(bool inDegrees = true) const;

  // Make 2 vectors into 3 orthogonal vectors
  static std::vector<V3D> makeVectorsOrthogonal(std::vector<V3D> &vectors);

  // Send to a stream
  void printSelf(std::ostream &) const;
  void readPrinted(std::istream &);
  void read(std::istream &);
  void write(std::ostream &) const;
  std::string toString() const;
  void fromString(const std::string &str);

  double volume() const {
    return fabs(m_vector.prod());
  } ///< Calculate the volume of a cube X*Y*Z

  int reBase(const V3D &, const V3D &,
             const V3D &); ///<rebase to new basis vector
  int masterDir(const double Tol =
                    1e-3) const; ///< Determine if there is a master direction
  bool nullVector(const double Tol = 1e-3) const {
    return m_vector.isZero(Tol);
  } ///< Determine if the point is null
  bool coLinear(const V3D &, const V3D &) const;

  void saveNexus(::NeXus::File *file, const std::string &name) const;
  void loadNexus(::NeXus::File *file, const std::string &name);

  template <typename T>
  friend V3D operator*(const Eigen::Matrix<T, 3, 3> &matrix, const V3D &vector);

private:
  Eigen::Vector3d m_vector;
};

// Overload operator <<
MANTID_KERNEL_DLL std::ostream &operator<<(std::ostream &, const V3D &);
MANTID_KERNEL_DLL std::istream &operator>>(std::istream &, V3D &);

template <typename T>
V3D operator*(const Eigen::Matrix<T, 3, 3> &matrix, const V3D &vector) {
  return V3D(matrix.template cast<double>() * vector.m_vector);
}

} // Namespace Kernel
} // Namespace Mantid

#endif /*MANTID_KERNEL_V3D_H_*/
