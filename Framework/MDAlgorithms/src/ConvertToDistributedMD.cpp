#include <MantidAPI/IEventWorkspace.h>
#include "MantidMDAlgorithms/ConvertToDistributedMD.h"
#include "MantidMDAlgorithms/EventToMDEventConverter.h"
#include "MantidKernel/PropertyWithValue.h"
#include "MantidKernel/make_unique.h"



namespace Mantid {
namespace MDAlgorithms {

using Mantid::Kernel::Direction;
using Mantid::API::WorkspaceProperty;
using namespace Mantid::Kernel;


// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(ConvertToDistributedMD)

//----------------------------------------------------------------------------------------------

/// Algorithms name for identification. @see Algorithm::name
const std::string ConvertToDistributedMD::name() const { return "ConvertToDistributedMD"; }

/// Algorithm's version for identification. @see Algorithm::version
int ConvertToDistributedMD::version() const { return 1; }

/// Algorithm's category for identification. @see Algorithm::category
const std::string ConvertToDistributedMD::category() const {
  return "TODO: FILL IN A CATEGORY";
}

/// Algorithm's summary for use in the GUI and help. @see Algorithm::summary
const std::string ConvertToDistributedMD::summary() const {
  return "TODO: FILL IN A SUMMARY";
}

//----------------------------------------------------------------------------------------------
/** Initialize the algorithm's properties.
 */
void ConvertToDistributedMD::init() {
  declareProperty(
      Kernel::make_unique<WorkspaceProperty<API::IEventWorkspace>>("InputWorkspace", "",
                                                             Direction::Input),
      "An input workspace.");
  declareProperty(
      Kernel::make_unique<WorkspaceProperty<API::Workspace>>("OutputWorkspace", "",
                                                             Direction::Output),
      "An output workspace.");

  declareProperty(make_unique<Mantid::Kernel::PropertyWithValue<bool>>("LorentzCorrection",
                      false, Direction::Input),
                  "Correct the weights of events by multiplying by the Lorentz "
                      "formula: sin(theta)^2 / lambda^4");

  declareProperty(
      make_unique<PropertyWithValue<double>>("Fraction", 0.01f),
      "Fraction of pulse time that should be used to build the initial data set.\n");

}

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void ConvertToDistributedMD::exec() {
  // Get the users's inputs
  Mantid::DataObjects::EventWorkspace_sptr inputWorkspace = getProperty("InputWorkspace");
  double fraction = getProperty("Fraction");

  // Get a n-percent fraction
  EventToMDEventConverter converter;
  auto nPercentEvents = converter.getEvents(*inputWorkspace, fraction, QFrame::QLab);

  // Build send all events to the master rank

  // Build the box structure on the master rank

  // Serialize the box structure

  // Share the box structure with all other ranks

  // Deserialize the box structure

  //


}

} // namespace MDAlgorithms
} // namespace Mantid
