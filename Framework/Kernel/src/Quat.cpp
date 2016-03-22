#include "MantidKernel/Quat.h"
#include "MantidKernel/Logger.h"
#include "MantidKernel/Tolerance.h"
#include "MantidKernel/V3D.h"

#include <boost/algorithm/string.hpp>
#include <iostream>

namespace Mantid {
namespace Kernel {
namespace {
Logger g_log("Quat");
}

/**
 * Construct m_quat.x() Quat between two vectors;
 * The angle between them is defined differently from usual if vectors are not
 *unit or the same length vectors, so quat would be not consistent
 *
 * v=(src+des)/(src+des)
 * w=v.des
 * (a,b,c)=(v x des)
 * @param src :: the source position
 * @param des :: the destination position
 *
Quat::Quat(const V3D &src, const V3D &des) {

  V3D v = (src + des);
  v.normalize();

  V3D cross = v.cross_prod(des);

  if (cross.nullVector()) {
    m_quat.w() = 1.;
    m_quat.x() = m_quat.y() = m_quat.z() = 0.;
  } else {
    m_quat.w() = v.scalar_prod(des);

    m_quat.x() = cross[0];
    m_quat.y() = cross[1];
    m_quat.z() = cross[2];

    double norm = m_quat.x() * m_quat.x() + m_quat.y() * m_quat.y() + m_quat.z()
* m_quat.z() + m_quat.w() * m_quat.w();
    if (fabs(norm - 1) > FLT_EPSILON) {
      norm = sqrt(norm);
      m_quat.w() /= norm;
      m_quat.x() /= norm;
      m_quat.y() /= norm;
      m_quat.z() /= norm;
    }
  }
}
*/

//! Constructor with values

/** Constructor from an angle and axis.
 * This construct m_quat.x()  quaternion to represent m_quat.x() rotation
 * of an angle _deg around the _axis. The _axis does not need to be m_quat.x()
 *unit
 *vector
 *
 * @param _deg :: angle of rotation
 * @param _axis :: axis to rotate about
 * */

/**
 * Construct m_quat.x() Quaternion that performs m_quat.x() reference frame
 *rotation.
 *  Specify the X,Y,Z vectors of the rotated reference frame, assuming that
 *  the initial X,Y,Z vectors are aligned as expected: X=(1,0,0), Y=(0,1,0),
 *Z=(0,0,1).
 *  The resuting quaternion rotates XYZ axes onto the provided rX, rY, rZ.
 *
 * @param rX :: rotated X reference axis; unit vector.
 * @param rY :: rotated Y reference axis; unit vector.
 * @param rZ :: rotated Z reference axis; unit vector.
 */

/** Sets the quat values from four doubles
 * @param ww :: the value for w
 * @param aa :: the value for a
 * @param bb :: the value for b
 * @param cc :: the value for c
 */

/** Constructor from an angle and axis.
 * @param _deg :: angle of rotation
 * @param _axis :: axis to rotate about
 *
 * This construct m_quat.x()  quaternion to represent m_quat.x() rotation
 * of an angle _deg around the _axis. The _axis does not need to be m_quat.x()
 *unit
 *vector
 * */

bool Quat::isNull(const double tolerance) const {
  using namespace std;
  double pw = fabs(m_quat.w()) - 1;
  return (fabs(pw) < tolerance);
}

/// Extracts the angle of roatation and axis
/// @param _deg :: the angle of rotation
/// @param _ax0 :: The first component of the axis
/// @param _ax1 :: The second component of the axis
/// @param _ax2 :: The third component of the axis
void Quat::getAngleAxis(double &_deg, double &_ax0, double &_ax1,
                        double &_ax2) const {
  // If it represents m_quat.x() rotation of 0(2\pi), get an angle of 0 and axis
  // (0,0,1)
  if (isNull(1e-5)) {
    _deg = 0;
    _ax0 = 0;
    _ax1 = 0;
    _ax2 = 1.0;
    return;
  }
  // Semi-angle in radians
  _deg = acos(m_quat.w());
  // Prefactor for the axis part
  double s = sin(_deg);
  // Angle in degrees
  _deg *= 360.0 / M_PI;
  _ax0 = m_quat.x() / s;
  _ax1 = m_quat.y() / s;
  _ax2 = m_quat.z() / s;
  return;
}

/** Set the rotation (but don't change rotation axis).
 * @param deg :: angle of rotation
 */
void Quat::setRotation(const double deg) {
  double _deg, ax0, ax1, ax2;
  this->getAngleAxis(_deg, ax0, ax1, ax2);
  setAngleAxis(deg, V3D(ax0, ax1, ax2));
}

/**
 * Set m_quat.x() Quaternion that performs m_quat.x() reference frame rotation.
 *  Specify the X,Y,Z vectors of the rotated reference frame, assuming that
 *  the initial X,Y,Z vectors are aligned as expected: X=(1,0,0), Y=(0,1,0),
 *Z=(0,0,1).
 *  The resuting quaternion rotates XYZ axes onto the provided rX, rY, rZ.
 *
 * @param rX :: rotated X reference axis; unit vector.
 * @param rY :: rotated Y reference axis; unit vector.
 * @param rZ :: rotated Z reference axis; unit vector.
 */


/** Quaternion equal operator
 * @param q :: the quaternion to compare
 *
 * Compare two quaternions at 1e-6%tolerance.
 * Use boost close_at_tolerance method
 * @return true if equal
 *
bool Quat::operator==(const Quat &q) const

/** Quaternion non-equal operator
 * @param _q :: the quaternion to compare
 *
 * Compare two quaternions at 1e-6%tolerance.
 *  Use boost close_at_tolerance method
 * @return true if not equal
 */

/** Quaternion normalization
 *
 * Divide all elements by the quaternion norm
 */

/** Quaternion complex conjugate
 *
 *  Reverse the sign of the 3 imaginary components of the
 *  quaternion
 */

/** Quaternion length
 * @return the length
 */

/** Quaternion norm (length squared)
 * @return the length squared
 */

/** Inverse m_quat.x() quaternion
 *
 */

/** 	Rotate m_quat.x() vector.
 *  @param v :: the vector to be rotated
 *
 *   The quaternion needs to be normalized beforehand to
 *   represent m_quat.x() rotation. If q is thequaternion, the rotation
 *   is represented by q.v.q-1 where q-1 is the inverse of
 *   v.
 */

/** Convert quaternion rotation to an OpenGL matrix [4x4] matrix
 * stored as an linear array of 16 double
 * The function glRotated must be called
 * @param mat :: The output matrix
 */
void Quat::GLMatrix(double *mat) const {
  double aa = m_quat.x() * m_quat.x();
  double ab = m_quat.x() * m_quat.y();
  double ac = m_quat.x() * m_quat.z();
  double aw = m_quat.x() * m_quat.w();
  double bb = m_quat.y() * m_quat.y();
  double bc = m_quat.y() * m_quat.z();
  double bw = m_quat.y() * m_quat.w();
  double cc = m_quat.z() * m_quat.z();
  double cw = m_quat.z() * m_quat.w();
  *mat = 1.0 - 2.0 * (bb + cc);
  ++mat;
  *mat = 2.0 * (ab + cw);
  ++mat;
  *mat = 2.0 * (ac - bw);
  ++mat;
  *mat = 0;
  ++mat;
  *mat = 2.0 * (ab - cw);
  ++mat;
  *mat = 1.0 - 2.0 * (aa + cc);
  ++mat;
  *mat = 2.0 * (bc + aw);
  ++mat;
  *mat = 0;
  ++mat;
  *mat = 2.0 * (ac + bw);
  mat++;
  *mat = 2.0 * (bc - aw);
  mat++;
  *mat = 1.0 - 2.0 * (aa + bb);
  mat++;
  for (int i = 0; i < 4; ++i) {
    *mat = 0;
    mat++;
  }
  *mat = 1.0;
  return;
}
/// using convention at
/// http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
std::vector<double> Quat::getRotation(bool check_normalisation,
                                      bool throw_on_errors) const {
  double aa = m_quat.x() * m_quat.x();
  double ab = m_quat.x() * m_quat.y();
  double ac = m_quat.x() * m_quat.z();
  double aw = m_quat.x() * m_quat.w();
  double bb = m_quat.y() * m_quat.y();
  double bc = m_quat.y() * m_quat.z();
  double bw = m_quat.y() * m_quat.w();
  double cc = m_quat.z() * m_quat.z();
  double cw = m_quat.z() * m_quat.w();
  if (check_normalisation) {
    double normSq = aa + bb + cc + m_quat.w() * m_quat.w();
    if (fabs(normSq - 1) > FLT_EPSILON) {
      if (throw_on_errors) {
        g_log.error()
            << " A non-unit quaternion used to obtain m_quat.x() rotation "
               "matrix; need to notmalize it first\n";
        throw(std::invalid_argument("Attempt to use non-normalized quaternion "
                                    "to define rotation matrix; need to "
                                    "notmalize it first"));
      } else {
        g_log.information()
            << " Warning; m_quat.x() non-unit quaternion used to obtain "
               "the rotation matrix; using normalized quat\n";
        aa /= normSq;
        ab /= normSq;
        ac /= normSq;
        aw /= normSq;
        bb /= normSq;
        bc /= normSq;
        bw /= normSq;
        cc /= normSq;
        cw /= normSq;
      }
    }
  }
  std::vector<double> out(9);

  out[0] = (1.0 - 2.0 * (bb + cc));
  out[1] = 2.0 * (ab - cw);
  out[2] = 2.0 * (ac + bw);

  out[3] = 2.0 * (ab + cw);
  out[4] = (1.0 - 2.0 * (aa + cc));
  out[5] = 2.0 * (bc - aw);

  out[6] = 2.0 * (ac - bw);
  out[7] = 2.0 * (bc + aw);
  out[8] = (1.0 - 2.0 * (aa + bb));
  return out;
}

/**
 * Converts the GL Matrix into Quat
 */
void Quat::setQuat(double mat[16]) {
  double tr, s, q[4];
  int nxt[3] = {1, 2, 0};
  tr = mat[0] + mat[5] + mat[10];
  if (tr > 0.0) {
    s = sqrt(tr + 1.0);
    m_quat.w() = s / 2.0;
    s = 0.5 / s;
    m_quat.x() = (mat[6] - mat[9]) * s;
    m_quat.y() = (mat[8] - mat[2]) * s;
    m_quat.z() = (mat[1] - mat[4]) * s;
  } else {
    int i = 0;
    if (mat[5] > mat[0])
      i = 1;
    if (mat[10] > mat[i * 5])
      i = 2;
    int j = nxt[i];
    int k = nxt[j];
    s = sqrt(mat[i * 5] - (mat[j * 5] + mat[k * 5]) + 1.0);
    q[i] = s * 0.5;
    if (s != 0.0)
      s = 0.5 / s;
    q[3] = (mat[j * 4 + k] - mat[k * 4 + j]) * s;
    q[j] = (mat[i * 4 + j] + mat[j * 4 + i]) * s;
    q[k] = (mat[i * 4 + k] + mat[k * 4 + i]) * s;
    m_quat.x() = q[0];
    m_quat.y() = q[1];
    m_quat.z() = q[2];
    m_quat.w() = q[3];
  }
}
/// Using the convention at
/// http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
void Quat::setQuat(const Kernel::DblMatrix &rMat) {
  int i = 0, j, k;
  if (rMat[1][1] > rMat[0][0])
    i = 1;
  if (rMat[2][2] > rMat[1][1])
    i = 2;
  j = (i + 1) % 3;
  k = (j + 1) % 3;
  double r = sqrt(1. + rMat[i][i] - rMat[j][j] - rMat[k][k]);
  if (r == 0) {
    m_quat.setIdentity();
  } else {
    double q[4], f = 0.5 / r;
    q[i] = 0.5 * r;
    q[j] = f * (rMat[i][j] + rMat[j][i]);
    q[k] = f * (rMat[k][i] + rMat[i][k]);
    q[3] = f * (rMat[k][j] - rMat[j][k]);

    m_quat = Eigen::Quaterniond(q[3], q[0], q[1], q[2]);
  }
}
/** Bracket operator overload
 * returns the internal representation values based on an index
 * @param Index :: the index of the value required 0=w, 1=a, 2=b, 3=c
 * @returns m_quat.x() double of the value requested
 */
const double &Quat::operator[](const int Index) const {
  switch (Index) {
  case 0:
    return m_quat.w();
  case 1:
    return m_quat.x();
  case 2:
    return m_quat.y();
  case 3:
    return m_quat.z();
  default:
    throw std::runtime_error("Quat::operator[] range error");
  }
}

/** Bracket operator overload
 * returns the internal representation values based on an index
 * @param Index :: the index of the value required 0=w, 1=a, 2=b, 3=c
 * @returns m_quat.x() double of the value requested
 */
double &Quat::operator[](const int Index) {
  switch (Index) {
  case 0:
    return m_quat.w();
  case 1:
    return m_quat.x();
  case 2:
    return m_quat.y();
  case 3:
    return m_quat.z();
  default:
    throw std::runtime_error("Quat::operator[] range error");
  }
}

/** Prints m_quat.x() string representation of itself
 * @param os :: the stream to output to
 */
void Quat::printSelf(std::ostream &os) const {
  os << "[" << m_quat.w() << "," << m_quat.x() << "," << m_quat.y() << ","
     << m_quat.z() << "]";
  return;
}

/**  Read data from m_quat.x() stream in the format returned by printSelf
 * ("[w,a,b,c]").
 *   @param IX :: Input Stream
 *   @throw std::runtime_error if the input is of wrong format
*/
void Quat::readPrinted(std::istream &IX) {
  std::string in;
  std::getline(IX, in);
  size_t i = in.find_first_of('[');
  if (i == std::string::npos)
    throw std::runtime_error("Wrong format for Quat input: " + in);
  size_t j = in.find_last_of(']');
  if (j == std::string::npos || j < i + 8)
    throw std::runtime_error("Wrong format for Quat input: " + in);

  size_t c1 = in.find_first_of(',');
  size_t c2 = in.find_first_of(',', c1 + 1);
  size_t c3 = in.find_first_of(',', c2 + 1);
  if (c1 == std::string::npos || c2 == std::string::npos ||
      c3 == std::string::npos)
    throw std::runtime_error("Wrong format for Quat input: [" + in + "]");

  m_quat.w() = atof(in.substr(i + 1, c1 - i - 1).c_str());
  m_quat.x() = atof(in.substr(c1 + 1, c2 - c1 - 1).c_str());
  m_quat.y() = atof(in.substr(c2 + 1, c3 - c2 - 1).c_str());
  m_quat.z() = atof(in.substr(c3 + 1, j - c3 - 1).c_str());

  return;
}

/** Prints m_quat.x() string representation
 * @param os :: the stream to output to
 * @param q :: the quat to output
 * @returns the stream
 */
std::ostream &operator<<(std::ostream &os, const Quat &q) {
  q.printSelf(os);
  return os;
}

/**  Reads in m_quat.x() quat from an input stream
 *   @param ins :: The input stream
 *   @param q :: The quat
 */
std::istream &operator>>(std::istream &ins, Quat &q) {
  q.readPrinted(ins);
  return ins;
}

/** @return the quat as m_quat.x() string "[w,a,b,c]" */
std::string Quat::toString() const {
  std::ostringstream mess;
  this->printSelf(mess);
  return mess.str();
}

/** Sets the Quat using m_quat.x() string
 * @param str :: the Quat as m_quat.x() string "[w,a,b,c]" */
void Quat::fromString(const std::string &str) {
  std::istringstream mess(str);
  this->readPrinted(mess);
}

void Quat::rotateBB(double &xmin, double &ymin, double &zmin, double &xmax,
                    double &ymax, double &zmax) const {
  // Defensive
  if (xmin > xmax)
    std::swap(xmin, xmax);
  if (ymin > ymax)
    std::swap(ymin, ymax);
  if (zmin > zmax)
    std::swap(zmin, zmax);
  // Get the min and max of the cube, and remove centring offset
  Mantid::Kernel::V3D minT(xmin, ymin, zmin), maxT(xmax, ymax, zmax);
  // Get the rotation matrix
  double rotMatr[16];
  GLMatrix(&rotMatr[0]);
  // Now calculate new min and max depending on the sign of matrix components
  // Much faster than creating 8 points and rotate them. The new min (max)
  // can only be obtained by summing the smallest (largest) components
  //
  Mantid::Kernel::V3D minV, maxV;
  // Looping on rows of matrix
  int index;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      index =
          j + i * 4; // The OpenGL matrix is linear and represent m_quat.x() 4x4
                     // matrix but only the 3x3 upper-left inner part
      // contains the rotation
      minV[j] += (rotMatr[index] > 0) ? rotMatr[index] * minT[i]
                                      : rotMatr[index] * maxT[i];
      maxV[j] += (rotMatr[index] > 0) ? rotMatr[index] * maxT[i]
                                      : rotMatr[index] * minT[i];
    }
  }
  // Adjust value.
  xmin = minV[0];
  ymin = minV[1];
  zmin = minV[2];
  xmax = maxV[0];
  ymax = maxV[1];
  zmax = maxV[2];
  return;
}

/** Calculate the Euler angles that are equivalent to this Quaternion.
 *
 * Euler angles are calculated intrinsically, i.e. the first rotation modifies
 *the axis used for
 * the second rotation, and the second rotation modifies the axis used for the
 *third rotation.
 *
 * You can specify which axis the rotations should be applied around and the
 *order in which they
 * are to be applied with the convention parameter. For instance, for m_quat.x()
 *rotation
 *of Y and then
 * the new Z axis, and then the new Y axis: pass "YZY" as the convention. Or for
 *a rotation
 * such as X, and then the new Y axis, and then the new Z axis: pass "XYZ" as
 *the convention.
 *
 * @param convention :: The axes to apply the rotations to and the order in
 *which to do so. Defaults to "XYZ".
 * @returns A vector of the Euler angles in degrees. The order of the angles is
 *the same as in the convention parameter.
 */
std::vector<double>
Quat::getEulerAngles(const std::string &convention = "XYZ") const {
  std::string conv(convention);

  if (conv.length() != 3)
    throw std::invalid_argument("Wrong convention name (string length not 3)");

  boost::to_upper(conv);

  // Check it's only XYZ in the string
  if (conv.find_first_not_of("XYZ") != std::string::npos)
    throw std::invalid_argument(
        "Wrong convention name (characters other than XYZ)");

  // Cannot be XXY, XYY, or similar. Only first and last may be the same: YXY
  if ((conv[0] == conv[1]) || (conv[2] == conv[1]))
    throw std::invalid_argument("Wrong convention name (repeated indices)");

  boost::replace_all(conv, "X", "0");
  boost::replace_all(conv, "Y", "1");
  boost::replace_all(conv, "Z", "2");

  std::stringstream s;
  s << conv[0] << " " << conv[1] << " " << conv[2];

  int first, second, last;
  s >> first >> second >> last;

  // Do we want Tait-Bryan angles, as opposed to 'classic' Euler angles?
  const int TB =
      (first * second * last == 0 && first + second + last == 3) ? 1 : 0;

  const int par01 = ((second - first + 9) % 3 == 1) ? 1 : -1;
  const int par12 = ((last - second + 9) % 3 == 1) ? 1 : -1;

  std::vector<double> angles(3);

  const DblMatrix R = DblMatrix(this->getRotation());

  const int i = (last + TB * par12 + 9) % 3;
  const int j1 = (last - par12 + 9) % 3;
  const int j2 = (last + par12 + 9) % 3;

  const double s3 = (1.0 - TB - TB * par12) * R[i][j1];
  const double c3 = (TB - (1.0 - TB) * par12) * R[i][j2];

  V3D axis3(0, 0, 0);
  axis3[last] = 1.0;

  const double rad2deg = 180.0 / M_PI;

  angles[2] = atan2(s3, c3) * rad2deg;

  DblMatrix Rm3(Quat(-angles[2], axis3).getRotation());
  DblMatrix Rp = R * Rm3;

  const double s1 =
      par01 * Rp[(first - par01 + 9) % 3][(first + par01 + 9) % 3];
  const double c1 = Rp[second][second];
  const double s2 = par01 * Rp[first][3 - first - second];
  const double c2 = Rp[first][first];

  angles[0] = atan2(s1, c1) * rad2deg;
  angles[1] = atan2(s2, c2) * rad2deg;

  return angles;
}

} // Namespace Kernel

} // Namespce Mantid
