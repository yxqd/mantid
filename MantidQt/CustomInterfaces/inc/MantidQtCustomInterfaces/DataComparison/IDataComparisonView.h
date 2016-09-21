#ifndef MANTIDQTCUSTOMINTERFACES_IDATACOMPARISONVIEW_H_
#define MANTIDQTCUSTOMINTERFACES_IDATACOMPARISONVIEW_H_

#include "MantidQtCustomInterfaces/DllConfig.h"
#include <string>

namespace MantidQt {
namespace CustomInterfaces {

/**
Interface for the view as in the MVP pattern. Methods declared here will be used
by the presenter to act on the view.

Copyright &copy; 2016 ISIS Rutherford Appleton Laboratory, NScD
Oak Ridge National Laboratory & European Spallation Source

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

class DLLExport IDataComparisonView {
public:
  /// Constructor
  IDataComparisonView(){};
  /// Destructor
  virtual ~IDataComparisonView(){};

  /// Print error message
  virtual void printError(const std::string &message) = 0;
  /// Print information message
  virtual void printInformation(const std::string &message) = 0;
  /// Pring debug message
  virtual void printDebug(const std::string &message) = 0;
};
}
}

#endif /* MANTIDQTCUSTOMINTERFACES_IDATACOMPARISONVIEW_H_ */
