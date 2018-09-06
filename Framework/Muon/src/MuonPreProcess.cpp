#include "MantidMuon/MuonPreProcess.h"
#include "MantidAPI/Algorithm.h"
#include "MantidAPI/AlgorithmManager.h"
#include "MantidDataObjects/TableWorkspace.h"

using namespace Mantid::API;
using namespace Mantid::DataObjects;
using namespace Mantid::Kernel;

namespace Mantid {
namespace Muon {

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(MuonPreProcess)

void MuonPreProcess::init() {
  std::string emptyString("");

  declareProperty(
      Mantid::Kernel::make_unique<WorkspaceProperty<Workspace>>(
          "InputWorkspace", "", Direction::Input, PropertyMode::Mandatory),
      "Input workspace containing data from detectors that the "
      "grouping/pairing will be applied to.");

  declareProperty("TimeMin", 0.0, "Start time for the data in micro seconds.",
                  Direction::Input);

  declareProperty("TimeMax", 32.0, "End time for the data in micro seconds.",
                  Direction::Input);

  declareProperty("RebinArgs", emptyString,
                  "Rebin arguments. No rebinning if left empty.",
                  Direction::Input);

  declareProperty("TimeOffset", 0.0,
                  "Shift the times of all data by a fixed amount (in micro "
                  "seconds). The value given corresponds to the bin that will "
                  "become 0.0 seconds.",
                  Direction::Input);

  declareProperty(
      make_unique<WorkspaceProperty<TableWorkspace>>(
          "DeadTimeTable", "", Direction::Input, PropertyMode::Optional),
      "Table with dead time information, used to apply dead time correction.");

  // Perform Group Associations.

  std::string analysisGrp("Analysis Options");
  setPropertyGroup("TimeMin", analysisGrp);
  setPropertyGroup("TimeMax", analysisGrp);
  setPropertyGroup("RebinArgs", analysisGrp);
  setPropertyGroup("TimeOffset", analysisGrp);
  setPropertyGroup("DeadTimeTable", analysisGrp);
}

void MuonPreProcess::exec() {}

} // namespace Muon
} // namespace Mantid
