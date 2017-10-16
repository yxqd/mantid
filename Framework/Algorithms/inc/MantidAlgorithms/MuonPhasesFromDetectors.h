#ifndef MANTID_ALGORITHMS_MUONPHASESFROMDETECTORS_H_
#define MANTID_ALGORITHMS_MUONPHASESFROMDETECTORS_H_

#include "MantidAlgorithms/DllConfig.h"
#include "MantidAPI/Algorithm.h"
#include "MantidAPI/ITableWorkspace_fwd.h"
#include "MantidAPI/WorkspaceGroup_fwd.h"
#include "MantidGeometry/IDTypes.h"

namespace Mantid {
namespace Indexing {
class SpectrumNumber;
}
namespace Algorithms {

/** MuonPhasesFromDetectors : Calculates asymmetry and phase for each spectra in a
  workspace based on the detector positions

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
class DLLExport MuonPhasesFromDetectors : public API::Algorithm {
public:
  /// Algorithm's name for identification overriding a virtual method
  const std::string name() const override { return "MuonPhasesFromDetectors"; }
  /// Summary of algorithms purpose
  const std::string summary() const override {
    return "Calculates the asymmetry and phase for each detector in a "
           "workspace based on the detector positions.";
  }

  /// Algorithm's version for identification overriding a virtual method
  int version() const override { return 1; }
  /// Algorithm's category for identification overriding a virtual method
  const std::string category() const override { return "Muon"; }

protected:
  /// Validate the inputs
  //std::map<std::string, std::string> validateInputs() override;

private:
  /// Initialise the algorithm
  void init() override;
  /// Execute the algorithm
  void exec() override;

  /// Report progress in GUI
  void reportProgress(const int thisSpectrum, const int totalSpectra);
  /// Pointer to input workspace
  API::MatrixWorkspace_sptr m_inputWS;
};
} // namespace Algorithms
} // namespace Mantid

#endif /* MANTID_ALGORITHMS_CALMUONDETECTORPHASES_H_ */
