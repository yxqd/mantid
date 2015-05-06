#ifndef MANTID_SINQ_POLDISPECTRUMCALIBRATIONFUNCTION_H_
#define MANTID_SINQ_POLDISPECTRUMCALIBRATIONFUNCTION_H_

#include "MantidSINQ/DllConfig.h"
#include "MantidSINQ/PoldiUtilities/PoldiSpectrumDomainFunction.h"

namespace Mantid {
namespace Poldi {

/** PoldiSpectrumCalibrationFunction

      @author Michael Wedel, Paul Scherrer Institut - SINQ
      @date 01/05/2015

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
class MANTID_SINQ_DLL PoldiSpectrumCalibrationFunction
    : public PoldiSpectrumDomainFunction {
public:
  PoldiSpectrumCalibrationFunction();
  virtual ~PoldiSpectrumCalibrationFunction() {}

  virtual std::string name() const {
    return "PoldiSpectrumCalibrationFunction";
  }

protected:
  void init();

  void
  functionModificationPreHook(const Poldi2DHelper_sptr &poldi2DHelper) const;
  void
  functionModificationPostHook(const Poldi2DHelper_sptr &poldi2DHelper) const;

  void setPeakCenter(double newCenter, double chopperOffset) const;
};

} // namespace Poldi
} // namespace Mantid

#endif /* MANTID_SINQ_POLDISPECTRUMCALIBRATIONFUNCTION_H_ */
