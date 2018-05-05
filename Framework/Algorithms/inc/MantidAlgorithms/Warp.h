#ifndef MANTID_ALGORITHMS_WARP_H_
#define MANTID_ALGORITHMS_WARP_H_

#include "MantidAPI/Algorithm.h"
#include "MantidAlgorithms/DllConfig.h"

#include "MantidGeometry/muParser_Silent.h"

#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/ring.hpp>

namespace Mantid {
namespace Algorithms {

		/** Warps the horizontal and vertical coordinates of a histogram
		workspace's bins according to a given expression.

		Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
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

class MANTID_ALGORITHMS_DLL Warp : public API::Algorithm {
public:
	/// A 2D point type.
	using Point = boost::geometry::model::d2::point_xy<double>;
  /// Clockwise, non-closed ring type.
  using Ring = boost::geometry::model::ring<Point, true, false, std::vector>;
  /// Clockwise, non-closed polygon type.
  using Polygon = boost::geometry::model::polygon<Point, true, false, std::vector, std::vector>;
  /// A collection type of Rings.
  using PolygonVector = boost::geometry::model::multi_polygon<Polygon, std::vector>;
  /// A bounding box type.
  using Box = boost::geometry::model::box<Point>;
	Warp();
	const std::string name() const override;
	int version() const override;
	const std::string category() const override;
	const std::string summary() const override;
private:
	double m_t{0.};
	double m_x{0.};
	mu::Parser m_tExpressionParser{};
	mu::Parser m_xExpressionParser{};
	void init() override;
	void exec() override;
	double evaluate(mu::Parser &parser, const Point &p);
  Ring warp(const Ring &q);
};

} // namespace Algorithms
} // namespace Mantid
#endif // MANTID_ALGORITHMS_WARP_H_
