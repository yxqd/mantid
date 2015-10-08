#ifndef MANTID_DATAHANDLING_SAVENXTOMOFROMSTACKOFIMAGES_H_
#define MANTID_DATAHANDLING_SAVENXTOMOFROMSTACKOFIMAGES_H_

#include "MantidAPI/Algorithm.h"

namespace Mantid {
namespace DataHandling {

/**
  SaveNXTomoFromStackOfImages

  Copyright &copy; 2015 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
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
class DLLExport SaveNXTomoFromStackOfImages : public API::Algorithm {
public:
  /// Default constructor
  SaveNXTomoFromStackOfImages();
  /// Destructor
  virtual ~SaveNXTomoFromStackOfImages();

  /// Algorithm's name for identification overriding a virtual method
  virtual const std::string name() const;
  /// Algorithm's version for identification overriding a virtual method
  virtual int version() const;
  /// Algorithm's category for identification overriding a virtual method
  virtual const std::string category() const;

  /// Summary of algorithms purpose
  virtual const std::string summary() const;

private:
  /// Initialise the algorithm
  void init();

  /// Execute the algorithm
  void exec();
};

} // namespace DataHandling
} // namespace Mantid

#endif /* MANTID_DATAHANDLING_SAVENXTOMOFROMSTACKOFIMAGES_H_ */
