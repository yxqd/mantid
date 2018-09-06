#include "MantidMuon/MuonGroupingCounts.h"
#include "MantidAPI/Algorithm.h"
#include "MantidAPI/AlgorithmManager.h"
#include "MantidDataObjects/TableWorkspace.h"

using namespace Mantid::API;
using namespace Mantid::DataObjects;
using namespace Mantid::Kernel;

namespace Mantid {
namespace Muon {

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(MuonGroupingCounts)

void MuonGroupingCounts::init() {
  std::string emptyString("");

  declareProperty(Mantid::Kernel::make_unique<WorkspaceProperty<Workspace>>(
                      "InputWorkspace", emptyString, Direction::Input,
                      PropertyMode::Mandatory),
                  "Input workspace containing data from detectors which are to "
                  "be grouped.");

  declareProperty(Mantid::Kernel::make_unique<WorkspaceProperty<Workspace>>(
                      "OutputWorkspace", emptyString, Direction::Output),
                  "Output workspace which will hold the grouped data.");

  declareProperty("GroupName", emptyString,
                  "The name of the group. Must "
                  "contain at least one alphanumeric "
                  "character.",
                  Direction::Input);
  declareProperty("Grouping", std::to_string(1),
                  "The grouping of detectors, comma separated list of detector "
                  "IDs or hyphenated ranges of IDs.",
                  Direction::Input);

  declareProperty("SummedPeriods", std::to_string(1),
                  "A list of periods to sum in multiperiod data.",
                  Direction::Input);
  declareProperty("SubtractedPeriods", emptyString,
                  "A list of periods to subtract in multiperiod data.",
                  Direction::Input);

  // Perform Group Associations.

  std::string groupingGrp("Grouping Information");
  setPropertyGroup("GroupName", groupingGrp);
  setPropertyGroup("Grouping", groupingGrp);

  std::string periodGrp("Multi-period Data");
  setPropertyGroup("SummedPeriods", periodGrp);
  setPropertyGroup("SubtractedPeriods", periodGrp);
}

void MuonGroupingCounts::exec() {}

} // namespace Muon
} // namespace Mantid
