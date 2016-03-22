#ifndef MANTID_KERNEL_QUAT_H_
#define MANTID_KERNEL_QUAT_H_

#include "MantidKernel/DllConfig.h"
#include "MantidKernel/Matrix.h"
#include "MantidKernel/V3D.h"

#include <Eigen/Geometry>

namespace Mantid {
namespace Kernel {
// Forward declarations

/** @class Quat Quat.h Geometry/Quat.h
@brief Class for quaternions
@version 1.0
@author Laurent C Chapon, ISIS RAL
@date 10/10/2007

Templated class for quaternions.
Quaternions are the 3D generalization of complex numbers
Quaternions are used for roations in 3D spaces and
often implemented for computer graphics applications.
Quaternion can be written q=W+ai+bj+ck where
w is the scalar part, and a, b, c the 3 imaginary parts.
Quaternion multiplication is non-commutative.<br>
i*j=-j*i=k<br>
j*k=-k*j=i<br>
k*i=-i*k=j<br>
Rotation of an angle theta around a normalized axis (u,v,w) can be simply
written W=cos(theta/2), a=u*sin(theta/2), b=v*sin(theta/2), c=w*sin(theta/2)
This class support all arithmetic operations for quaternions

Copyright &copy; 2007 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
National Laboratory & European Spallation Source

This file is part of Mantid.

Mantid is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your ption) any later version.

Mantid is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

File change history is stored at: <https://github.com/mantidproject/mantid>
*/
class MANTID_KERNEL_DLL Quat {

public:
  inline Quat() : m_quat(Eigen::Quaterniond::Identity()) {}
  // direct quat definition
  inline Quat(const double _w, const double _a, const double _b,
              const double _c)
      : m_quat(Eigen::Quaterniond(_w, _a, _b, _c)) {}

  // * Construct a Quat between two vectors;
  // * The angle between them is defined differently from usual if vectors are
  // not unit or the same length vectors, so quat would be not consistent
  inline Quat(const V3D &src, const V3D &des)
      : m_quat(Eigen::Quaterniond::FromTwoVectors(src.getVector().normalized(),
                                                  des.getVector().normalized())
                   .normalized()) {}

  inline Quat(const V3D &rX, const V3D &rY, const V3D &rZ) {
    // Call the operator to do the setting
    this->operator()(rX, rY, rZ);
  }

  //! Set quaternion form an angle in degrees and an axis
  inline Quat(const double _deg, const V3D &_axis)
      : m_quat(Eigen::AngleAxisd(_deg / 180.0 * M_PI,
                                 _axis.getVector().normalized())) {}

  explicit Quat(const Eigen::Quaterniond &quat) : m_quat(quat) {}

  // set a quaternion from a rotational matrix;
  Quat(const DblMatrix &RotMat) { this->setQuat(RotMat); }
  inline void operator()(const Quat &q) { m_quat = q.m_quat; }
  inline void operator()(const double ww, const double aa, const double bb,
                         const double cc) {
    set(ww, aa, bb, cc);
  }
  void operator()(const double angle, const V3D &v) {
    m_quat = Eigen::Quaterniond(
        Eigen::AngleAxisd(angle / 180.0 * M_PI, v.getVector().normalized()));
  }

  void operator()(const V3D &rX, const V3D &rY, const V3D &rZ) {
    // The quaternion will combine two quaternions.
    UNUSED_ARG(rZ); // Avoid compiler warning

    Eigen::Quaterniond Q1 = Eigen::Quaterniond::FromTwoVectors(
        Eigen::Vector3d(1, 0, 0), rX.getVector());

    Eigen::Quaterniond Q2 = Eigen::Quaterniond::FromTwoVectors(
        Q1._transformVector(Eigen::Vector3d(0, 1, 0)), rY.getVector());

    m_quat = Q2 * Q1;
  }

  void set(const double ww, const double aa, const double bb, const double cc) {
    m_quat = Eigen::Quaterniond(ww, aa, bb, cc);
  }
  void setAngleAxis(const double _deg, const V3D &_axis) {
    m_quat = Eigen::Quaterniond(
        Eigen::AngleAxisd(_deg * M_PI / 180.0, _axis.getVector().normalized()));
  }
  void getAngleAxis(double &_deg, double &_ax0, double &_ax1,
                    double &ax2) const;
  std::vector<double> getEulerAngles(const std::string &convention) const;
  /// Set the rotation (both don't change rotation axis)
  void setRotation(const double deg);
  //! Norm of a quaternion
  double len() const { return m_quat.norm(); }
  //! Norm squared
  double len2() const { return m_quat.squaredNorm(); }
  //! Re-initialize to identity
  void init() { m_quat.setIdentity(); }
  //! Normalize
  void normalize() { m_quat.normalize(); }
  //! Take the complex conjugate
  void conjugate() { m_quat = m_quat.conjugate(); }
  //! Inverse a quaternion (in the sense of rotation inversion)
  void inverse() { m_quat = m_quat.inverse(); }
  //! Is the quaternion representing a null rotation
  bool isNull(const double tolerance = 0.001) const;
  //! Convert quaternion rotation to an OpenGL matrix [4x4] matrix
  //! stored as an linear array of 16 double
  //! The function glRotated must be called
  void GLMatrix(double *mat) const;
  //! returns the rotation matrix defined by this quaternion as an 9-point
  // vector representing M33 matrix
  //! (m33 is not used at the moment), if check_normalisation selected, verify
  // if the mod(quat) is indeed == 1 and throws otherwise.
  std::vector<double> getRotation(bool check_normalisation = false,
                                  bool throw_on_errors = false) const;
  //! Convert GL Matrix into Quat
  void setQuat(double mat[16]);
  //! Convert usual 3D rotation matrix into quat; Will throw if matirix is not
  // rotational;
  void setQuat(const DblMatrix &rMat);
  //! Rotate a vector
  void rotate(V3D &v) const { v = V3D(m_quat._transformVector(v.getVector())); }

  //! Taking two points defining a cuboid bounding box (xmin,ymin,zmin) and
  //(xmax,ymax,zmax)
  // which means implicitly that the cube edges are parallel to the axes,
  // find the smallest bounding box with the edges also parallel to the axes
  // after rotation of the object.
  void rotateBB(double &xmin, double &ymin, double &zmin, double &xmax,
                double &ymax, double &zmax) const;
  //! Overload operators
  inline Quat operator+(const Quat &_q) const {
    Quat out(*this);
    out += _q;

    return out;
  }

  inline Quat &operator+=(const Quat &_q) {
    m_quat.w() += _q.m_quat.w();
    m_quat.x() += _q.m_quat.x();
    m_quat.y() += _q.m_quat.y();
    m_quat.z() += _q.m_quat.z();
    return *this;
  }

  inline Quat operator-(const Quat &_q) const {
    Quat out(*this);
    out -= _q;

    return out;
  }

  inline Quat &operator-=(const Quat &_q) {
    m_quat.w() -= _q.m_quat.w();
    m_quat.x() -= _q.m_quat.x();
    m_quat.y() -= _q.m_quat.y();
    m_quat.z() -= _q.m_quat.z();
    return *this;
  }

  inline Quat operator*(const Quat &_q) const {
    Quat out(*this);
    out *= _q;
    return out;
  }

  inline Quat &operator*=(const Quat &_q) {
    m_quat *= _q.m_quat;
    return *this;
  }

  inline bool operator==(const Quat &q) const {
    using namespace std;
    return !(fabs(m_quat.w() - q.m_quat.w()) > Tolerance ||
             fabs(m_quat.x() - q.m_quat.x()) > Tolerance ||
             fabs(m_quat.y() - q.m_quat.y()) > Tolerance ||
             fabs(m_quat.z() - q.m_quat.z()) > Tolerance);
  }

  inline bool operator!=(const Quat &q) const { return !(this->operator==(q)); }
  const double &operator[](int) const;
  double &operator[](int);

  /** @name Element access. */
  //@{
  /// Access the real part
  inline double real() const { return m_quat.w(); }
  /// Access the coefficient of i
  inline double imagI() const { return m_quat.x(); }
  /// Access the coefficient of j
  inline double imagJ() const { return m_quat.y(); }
  /// Access the coefficient of k
  inline double imagK() const { return m_quat.z(); }
  //@}

  void printSelf(std::ostream &) const;
  void readPrinted(std::istream &);
  std::string toString() const;
  void fromString(const std::string &str);

private:
  Eigen::Quaterniond m_quat;
};

MANTID_KERNEL_DLL std::ostream &operator<<(std::ostream &, const Quat &);
MANTID_KERNEL_DLL std::istream &operator>>(std::istream &, Quat &q);

} // Namespace Mantid

} // Namespace Kernel

#endif /*MANTID_KERNEL_QUAT_H_*/
