#ifndef MANTID_GEOMETRY_BASEHKL_H_
#define MANTID_GEOMETRY_BASEHKL_H_

#include <Eigen/Core>

namespace Mantid {
namespace Geometry {

/** BaseHKL : TODO: DESCRIPTION

  Copyright &copy; 2016 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
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

  File change history is stored at: <https://github.com/mantidproject/mantid>
  Code Documentation is available at: <http://doxygen.mantidproject.org>
*/
template <typename NumericType, typename T> class BaseHKL {
public:
  BaseHKL(const NumericType &h, const NumericType &k, const NumericType l)
      : m_hklRepresentation(h, k, l) {}
  BaseHKL(const Eigen::Matrix<NumericType, 3, 1> &hkl)
      : m_hklRepresentation(hkl) {}

  const NumericType &h() const { return m_hklRepresentation(0); }
  const NumericType &k() const { return m_hklRepresentation(1); }
  const NumericType &l() const { return m_hklRepresentation(2); }

  const Eigen::Matrix<NumericType, 3, 1> &getRepresentation() const {
    return m_hklRepresentation;
  }

  T operator+(const T &other) const {
    return T(m_hklRepresentation + other.m_hklRepresentation);
  }

  T operator-(const T &other) const {
    return T(m_hklRepresentation - other.m_hklRepresentation);
  }

  bool operator==(const T &other) const {
    return (m_hklRepresentation - other.m_hklRepresentation).isZero();
  }

  bool operator!=(const T &other) const { return !this->operator==(other); }

private:
  Eigen::Matrix<NumericType, 3, 1> m_hklRepresentation;
};

class ProHKL : public BaseHKL<double, ProHKL> {
public:
  ProHKL(double h, double k, double l) : BaseHKL(h, k, l) {}
  ProHKL(const Eigen::Vector3d &hkl) : BaseHKL(hkl) {}
};

class FractionalHKL : public BaseHKL<double, FractionalHKL> {
public:
  FractionalHKL(double h, double k, double l) : BaseHKL(h, k, l) {}
  FractionalHKL(const Eigen::Vector3d &hkl) : BaseHKL(hkl) {}

  explicit FractionalHKL(const ProHKL &proHKL)
      : BaseHKL(proHKL.getRepresentation()) {}
};

class IntegerHKL : public BaseHKL<int, IntegerHKL> {
public:
  IntegerHKL(int h, int k, int l) : BaseHKL(h, k, l) {}
  IntegerHKL(const Eigen::Vector3i &hkl) : BaseHKL(hkl) {}

  explicit IntegerHKL(const ProHKL &proHKL)
      : BaseHKL(proHKL.getRepresentation().unaryExpr([](const double &el) {
          return std::round(el);
        }).cast<int>()) {}

  explicit operator ProHKL() const {
    return ProHKL(getRepresentation().cast<double>());
  }
};

} // namespace Geometry
} // namespace Mantid

#endif /* MANTID_GEOMETRY_BASEHKL_H_ */
