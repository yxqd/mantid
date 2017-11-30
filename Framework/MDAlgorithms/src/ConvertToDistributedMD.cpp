#include <MantidAPI/IEventWorkspace.h>
#include "MantidMDAlgorithms/ConvertToDistributedMD.h"
#include "MantidMDAlgorithms/EventToMDEventConverter.h"
#include "MantidParallel/Communicator.h"
#include "MantidParallel/Collectives.h"


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

  // Get the preliminary box structure
  auto* boxes = getPreliminaryBoxStructure(nPercentEvents);
}


ConvertToDistributedMD::BoxStructure* ConvertToDistributedMD::getPreliminaryBoxStructure(ConvertToDistributedMD::MDEventList& mdEvents) const {
  ConvertToDistributedMD::BoxStructure* boxStructure= nullptr;

  const auto& communicator = this->communicator();
  if (communicator.rank() == 0) {
    // 1.b Receive the data from all the other ranks
    auto allEvents = receiveMDEvents(const Mantid::Parallel::Communicator& communicator, mdEvents)

    // 2. Build the box structure on the master rank

    // 4. Serialize the box structure

    // 5.a Broadcast serialized box structure to all other ranks

  } else {
    // 1.a Send the event data to the master rank
    sendMDEventsToMaster(communicator, mdEvents);

    // 5.b Receive the box structure from master

    // 6. Deserialize on all ranks (except for master)
  }



  return boxStructure;
}

void ConvertToDistributedMD::sendMDEventsToMaster(const Mantid::Parallel::Communicator& communicator,
                                                  const MDEventList& mdEvents) const {
  // Send signal array and qx qy qz arrays separately
  // This could also be separate arrays

  // Send the totalNumberOfEvents
  auto totalNumberEvents = mdEvents.size();
  gather(communicator, totalNumberEvents, 0);

  // Build up the vectors we want to send, i.e. signal and position;
  // TODO: use boost::serialization
  auto numberOfDimensions = mdEvents.back().getNumDims();
  std::vector<float> signals(totalNumberEvents);
  std::vector<coord_t> positions(totalNumberEvents*numberOfDimensions);
  for (const auto& event : mdEvents) {
    signals.emplace_back(event.getSignal());
    auto position = event.get
  }

  // Send the position array


}


ConvertToDistributedMD::MDEventList ConvertToDistributedMD::receiveMDEventsOnMaster(const Mantid::Parallel::Communicator& communicator,
                                                                            const MDEventList& mdEvents) const {
  // Get the totalNumberOfEvents for each rank
  std::vector<size_t> numberOfEventsPerRank;
  auto masterNumberOfEvents = mdEvents.size();
  gather(communicator, masterNumberOfEvents, numberOfEventsPerRank, 0);

  // Receive the signal arrays. To do this we need to:
  // 1. Create a sufficiently large buffer
  // 2. Determine the strides that each rank requires
  // 3. Apply gatherv

  // 1. Create a receive buffer for the signals
  auto totalNumberOfEvents = std::accumulate(numberOfEventsPerRank.begin(), numberOfEventsPerRank.end(), 0ul);
  std::vector<int> signals(totalNumberOfEvents);

  // 2. Create strides
  std::vector<int> strides(numberOfEventsPerRank.size());
  std::partial_sum(numberOfEventsPerRank.begin(), numberOfEventsPerRank.end(), strides.begin());
  auto correctFirstValue = [&numberOfEventsPerRank](int& value) {value -= numberOfEventsPerRank[0];};
  std::for_each(strides.begin(), strides.end(), correctFirstValue);

  // 3 Receive the data from the other ranks.

}



  } // namespace MDAlgorithms
} // namespace Mantid
