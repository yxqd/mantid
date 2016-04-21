#ifndef MANTID_CURVEFITTING_POLDICALIBRATIONFUNCTION_H_
#define MANTID_CURVEFITTING_POLDICALIBRATIONFUNCTION_H_

#include "MantidCurveFitting/Functions/Gaussian.h"

namespace Mantid {
namespace CurveFitting {
namespace Functions {
/** PoldiCalibrationProfile

  Instrument parameters at POLDI are calibrated using an additional parameter
  that changes the slope of Bragg lines in the POLDI 2D data. If the
  instrument is calibrated properly, the parameter is 0.

  Since silicon powder is used as standard, peak profiles are always Gaussian,
  so the Gaussian function from this module is re-used. The additional
  parameter depends on 2theta, so it can only be used in conjunction with the
  PoldiSpectrumCalibrationFunction in the SINQ module. If used otherwise
  it will behave like Gaussian with an additional parameter that does nothing.

  The attribute DeltaTheta has to be set in order to get the correct shift for
  the corresponding 2theta-value. The ChopperOffset attribute is another shift
  that is added to the peak centre after correcting for the 2theta-dependent
  shift, so that the real peak centre is:

    c_real = c_0 * [1.0 + shiftFactor(DeltaTheta)] + ChopperOffset

  The centre/setCentre methods are unmodified, so they behave exactly like in
  the Gaussian case, modifying only the PeakCentre parameter.

    @author Michael Wedel, Paul Scherrer Institut - SINQ
    @date 04/05/2015

  Copyright © 2015 PSI-NXMM

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
class DLLExport PoldiCalibrationProfile : public Gaussian {
public:
  ~PoldiCalibrationProfile() {}

  std::string name() const { return "PoldiCalibrationProfile"; }

  void functionLocal(double *out, const double *xValues,
                     const size_t nData) const;
  void functionDerivLocal(API::Jacobian *out, const double *xValues,
                          const size_t nData);

protected:
  double getShiftFactor() const;
  double getConstantFactor() const;

  void init();
};

}
} // namespace CurveFitting
} // namespace Mantid

#endif /* MANTID_CURVEFITTING_POLDICALIBRATIONFUNCTION_H_ */
