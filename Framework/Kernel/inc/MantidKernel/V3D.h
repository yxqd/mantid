#ifndef MANTID_KERNEL_V3D_H_
#define MANTID_KERNEL_V3D_H_

#include <cmath>
#include <cfloat>
#include <complex>
#include <vector>
#include "MantidKernel/DllConfig.h"
#include "MantidKernel/Matrix.h"
#include "MantidKernel/Tolerance.h"
#include <nexus/NeXusFile.hpp>

#include <Eigen/Core>

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
  V3D();
  V3D(const double, const double, const double);
  explicit V3D(const Eigen::Vector3d &vector);

  /// Convenience method for sorting list of V3D objects based on magnitude
  static bool CompareMagnitude(const Kernel::V3D &v1, const Kernel::V3D &v2);

  // explicit conversion into vector
  operator std::vector<double>() const {
    return std::vector<double>{m_vector.data(),
                               m_vector.data() + m_vector.size()};
  }

  // Arithemetic operators overloaded
  inline V3D operator+(const V3D &v) const {
    V3D out(*this);
    out += v;
    return out;
  }

  inline V3D &operator+=(const V3D &v) {
    m_vector += v.m_vector;
    return *this;
  }

  inline V3D operator-(const V3D &v) const {
    V3D out(*this);
    out -= v;
    return out;
  }

  inline V3D &operator-=(const V3D &v) {
    m_vector -= v.m_vector;
    return *this;
  }

  inline V3D operator*(const V3D &v) const {
    V3D out(*this);
    out *= v;
    return out;
  }

  inline V3D operator/(const V3D &v) const {
    V3D out(*this);
    out /= v;
    return out;
  }

  inline V3D &operator*=(const V3D &v) {
    m_vector = m_vector.cwiseProduct(v.m_vector);
    return *this;
  }

  inline V3D &operator/=(const V3D &v) {
    m_vector = m_vector.cwiseQuotient(v.m_vector);
    return *this;
  }

  // Scale
  V3D operator*(const double D) const;
  V3D &operator*=(const double D);
  V3D operator/(const double D) const;
  V3D &operator/=(const double D);
  // Simple Comparison
  inline bool operator==(const V3D &v) const {
    return (m_vector - v.m_vector).isZero(Kernel::Tolerance);
  }
  bool operator!=(const V3D &) const;
  bool operator<(const V3D &) const;
  bool operator>(const V3D &rhs) const;

  // Access
  // Setting x, y and z values
  void operator()(const double xx, const double yy, const double zz);
  void spherical(const double &R, const double &theta, const double &phi);
  void spherical_rad(const double &R, const double &polar,
                     const double &azimuth);
  void azimuth_polar_SNS(const double &R, const double &azimuth,
                         const double &polar);
  void setX(const double xx);
  void setY(const double yy);
  void setZ(const double zz);

  const double &X() const { return m_vector(0); } ///< Get x
  const double &Y() const { return m_vector(1); } ///< Get y
  const double &Z() const { return m_vector(2); } ///< Get z

  const double &operator[](const size_t Index) const;
  double &operator[](const size_t Index);

  void getSpherical(double &R, double &theta, double &phi) const;

  //      void rotate(const V3D&,const V3D&,const double);
  void rotate(const Matrix<double> &);

  void round();

  /// Make a normalized vector (return norm value)
  double normalize(); // Vec3D::makeUnit
  double norm() const;
  double norm2() const;
  /// transform vector into form, used to describe directions in
  /// crystallogaphical coodinate system
  double toMillerIndexes(double eps = 1.e-3);
  /// Scalar product
  inline double scalar_prod(const V3D &v) const {
    return m_vector.dot(v.m_vector);
  }
  /// Cross product
  inline V3D cross_prod(const V3D &v) const {
    return V3D(m_vector(1) * v.m_vector(2) - m_vector(2) * v.m_vector(1),
               m_vector(2) * v.m_vector(0) - m_vector(0) * v.m_vector(2),
               m_vector(0) * v.m_vector(1) - m_vector(1) * v.m_vector(0));
  }
  /// Distance (R) between two points defined as vectors
  double distance(const V3D &v) const;
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
    return fabs(m_vector(0) * m_vector(1) * m_vector(2));
  } ///< Calculate the volume of a cube X*Y*Z

  int reBase(const V3D &, const V3D &,
             const V3D &); ///<rebase to new basis vector
  int masterDir(const double Tol =
                    1e-3) const; ///< Determine if there is a master direction
  bool
  nullVector(const double Tol = 1e-3) const; ///< Determine if the point is null
  bool coLinear(const V3D &, const V3D &) const;

  void saveNexus(::NeXus::File *file, const std::string &name) const;
  void loadNexus(::NeXus::File *file, const std::string &name);

private:
  Eigen::Vector3d m_vector;
};

// Overload operator <<
MANTID_KERNEL_DLL std::ostream &operator<<(std::ostream &, const V3D &);
MANTID_KERNEL_DLL std::istream &operator>>(std::istream &, V3D &);

} // Namespace Kernel
} // Namespace Mantid

#endif /*MANTID_KERNEL_V3D_H_*/
