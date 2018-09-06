#include "MantidMuon/MuonPairingAsymmetry.h"
#include "MantidAPI/Algorithm.h"
#include "MantidAPI/AlgorithmManager.h"
#include "MantidDataObjects/TableWorkspace.h"
#include "MantidKernel/EnabledWhenProperty.h"

using namespace Mantid::API;
using namespace Mantid::DataObjects;
using namespace Mantid::Kernel;

namespace Mantid {
namespace Muon {

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(MuonPairingAsymmetry)

void MuonPairingAsymmetry::init() {
  std::string emptyString("");

  declareProperty(
      Mantid::Kernel::make_unique<WorkspaceProperty<WorkspaceGroup>>(
          "OutputWorkspaceGroup", emptyString, Direction::Output),
      "The workspace which will hold the results of the asymmetry "
      "calculation.");

  declareProperty("PairName", emptyString,
                  "The name of the pair. Must "
                  "contain at least one alphanumeric "
                  "character.",
                  Direction::Input);

  declareProperty("Alpha", 1.0,
                  "Alpha parameter used in the asymmetry calculation.",
                  Direction::Input);

  declareProperty("SpecifyGroupsManually", false,
                  "Specify the pair of groups manually using the raw data and "
                  "various optional parameters.");

  // Select groups via workspaces

  declareProperty(
      Mantid::Kernel::make_unique<WorkspaceProperty<MatrixWorkspace>>(
          "InputWorkspace1", emptyString, Direction::Input,
          PropertyMode::Optional),
      "Input workspace containing data from grouped detectors.");

  declareProperty(
      Mantid::Kernel::make_unique<WorkspaceProperty<MatrixWorkspace>>(
          "InputWorkspace2", emptyString, Direction::Input,
          PropertyMode::Optional),
      "Input workspace containing data from grouped detectors.");

  setPropertySettings("InputWorkspace1",
                      make_unique<Kernel::EnabledWhenProperty>(
                          "SpecifyGroupsManually", Kernel::IS_EQUAL_TO, "0"));
  setPropertySettings("InputWorkspace2",
                      make_unique<Kernel::EnabledWhenProperty>(
                          "SpecifyGroupsManually", Kernel::IS_EQUAL_TO, "0"));

  // Specify groups manually

  declareProperty(Mantid::Kernel::make_unique<WorkspaceProperty<Workspace>>(
                      "InputWorkspace", emptyString, Direction::Input,
                      PropertyMode::Optional),
                  "Input workspace containing data from detectors which are to "
                  "be grouped.");
  setPropertySettings("InputWorkspace",
                      make_unique<Kernel::EnabledWhenProperty>(
                          "SpecifyGroupsManually", Kernel::IS_EQUAL_TO, "1"));

  declareProperty("Group1", std::to_string(1),
                  "The grouping of detectors, comma separated list of detector "
                  "IDs or hyphenated ranges of IDs.",
                  Direction::Input);
  declareProperty("Group2", std::to_string(1),
                  "The grouping of detectors, comma separated list of detector "
                  "IDs or hyphenated ranges of IDs.",
                  Direction::Input);
  setPropertySettings("Group1",
                      make_unique<Kernel::EnabledWhenProperty>(
                          "SpecifyGroupsManually", Kernel::IS_EQUAL_TO, "1"));
  setPropertySettings("Group2",
                      make_unique<Kernel::EnabledWhenProperty>(
                          "SpecifyGroupsManually", Kernel::IS_EQUAL_TO, "1"));

  declareProperty("SummedPeriods", std::to_string(1),
                  "A list of periods to sum in multiperiod data.",
                  Direction::Input);
  setPropertySettings("SummedPeriods",
                      make_unique<Kernel::EnabledWhenProperty>(
                          "SpecifyGroupsManually", Kernel::IS_EQUAL_TO, "1"));

  declareProperty("SubtractedPeriods", emptyString,
                  "A list of periods to subtract in multiperiod data.",
                  Direction::Input);
  setPropertySettings("SubtractedPeriods",
                      make_unique<Kernel::EnabledWhenProperty>(
                          "SpecifyGroupsManually", Kernel::IS_EQUAL_TO, "1"));

  // Perform Group Associations.

  std::string workspaceGrp("Specify Group Workspaces");
  setPropertyGroup("InputWorkspace1", workspaceGrp);
  setPropertyGroup("InputWorkspace2", workspaceGrp);

  std::string manualGroupGrp("Specify Detector ID Groups Manually");
  setPropertyGroup("InputWorkspace", manualGroupGrp);
  setPropertyGroup("Group1", manualGroupGrp);
  setPropertyGroup("Group2", manualGroupGrp);

  std::string periodGrp("Multi-period Data");
  setPropertyGroup("SummedPeriods", periodGrp);
  setPropertyGroup("SubtractedPeriods", periodGrp);
}
void MuonPairingAsymmetry::exec() {}

} // namespace Muon
} // namespace Mantid
