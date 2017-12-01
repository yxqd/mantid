#include <MantidAPI/IEventWorkspace.h>
#include "MantidMDAlgorithms/ConvertToDistributedMD.h"
#include "MantidMDAlgorithms/EventToMDEventConverter.h"
#include "MantidParallel/Communicator.h"
#include "MantidParallel/Collectives.h"

#include <boost/mpi/request.hpp>

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
      Kernel::make_unique<WorkspaceProperty<DataObjects::EventWorkspace>>("InputWorkspace", "",
                                                             Direction::Input),
      "An input workspace.");
//  declareProperty(
//      Kernel::make_unique<WorkspaceProperty<DataObjects::EventWorkspace>("OutputWorkspace", "",
//                                                             Direction::Output),
//      "An output workspace.");

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
  // The generation of the distributed MDEventWorkspace (here only distributed box structure) has the following parts:
  // 1. Convert the data which corresponds to the first n percent of the total measurement
  // 2. Build a preliminary box structure on the master rank
  // 3. Determine how the data will be split based on the preliminary box structure and share this information with all
  //    ranks
  // 4. Share the preliminary box structure with all ranks
  // 5. Convert all events
  // 6. Disable the box controller and add events to the local box structure
  // 7. Send data from the each rank to the correct rank
  // 8. Enable the box controller and start splitting the data

  // Get the users's inputs
  Mantid::DataObjects::EventWorkspace_sptr inputWorkspace = getProperty("InputWorkspace");
  double fraction = getProperty("Fraction");

  // ----------------------------------------------------------
  // 1. Get a n-percent fraction
  // ----------------------------------------------------------
  EventToMDEventConverter converter;
  auto nPercentEvents = converter.getEvents(*inputWorkspace, fraction, QFrame::QLab);

  // -----------------------------------------------------------------
  // 2. + 3. = 4.  Get the preliminary box structure and the partition behaviour
  // -----------------------------------------------------------------
  auto* boxes = getPreliminaryBoxStructure(nPercentEvents);

  // ----------------------------------------------------------
  // 5. Convert all events
  // ----------------------------------------------------------

  // ----------------------------------------------------------
  // 6. Add the local data to the preliminary box structure
  // ----------------------------------------------------------

  // ----------------------------------------------------------
  // 7. Redistribute data
  // ----------------------------------------------------------

  // ----------------------------------------------------------
  // 8. Continue to split locally
  // ----------------------------------------------------------
}


ConvertToDistributedMD::BoxStructure* ConvertToDistributedMD::getPreliminaryBoxStructure(ConvertToDistributedMD::MDEventList& mdEvents) const {
  ConvertToDistributedMD::BoxStructure* boxStructure= nullptr;

  const auto& communicator = this->communicator();
  if (communicator.rank() == 0) {
    // 1.b Receive the data from all the other ranks
    auto allEvents = receiveMDEventsOnMaster(communicator, mdEvents);

    // 2. Build the box structure on the master rank
    //boxStructure = generatePreliminaryBoxStructure(allEvents);

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
  // Send the totalNumberOfEvents
  auto totalNumberEvents = mdEvents.size();
  gather(communicator, totalNumberEvents, 0);

  // Send the vector of md events
  communicator.send(0, 1, mdEvents.data(), mdEvents.size());
}


ConvertToDistributedMD::MDEventList ConvertToDistributedMD::receiveMDEventsOnMaster(const Mantid::Parallel::Communicator& communicator,
                                                                                    const MDEventList& mdEvents) const {
  // In order to get all relevant md events onto the master rank we need to:
  // 1. Inform the master rank about the number of events on all ranks (gather)
  // 2. Create a sufficiently large buff on master
  // 3. Determine the stride that each rank requires, ie which offset each rank requires
  // 4. Send the data to master. This would be a natural operation for gatherv, however this is not available in
  //    boost 1.58 (in 1.59+ it is).

  // 1. Get the totalNumberOfEvents for each rank
  std::vector<size_t> numberOfEventsPerRank;
  auto masterNumberOfEvents = mdEvents.size();
  gather(communicator, masterNumberOfEvents, numberOfEventsPerRank, 0);

  // 2. Create a buffer
  auto totalNumberOfEvents = std::accumulate(numberOfEventsPerRank.begin(), numberOfEventsPerRank.end(), 0ul);
  MDEventList totalEvents;
  totalEvents.reserve(totalNumberOfEvents);
  std::copy(mdEvents.begin(), mdEvents.end(), std::back_inserter(totalEvents));

  // 3. Determine the strides
  std::vector<int> strides;
  strides.reserve(numberOfEventsPerRank.size());
  strides.push_back(0);
  {
    std::vector<int> tempStrides(numberOfEventsPerRank.size()-1);
    std::partial_sum(numberOfEventsPerRank.begin(), numberOfEventsPerRank.end()-1, tempStrides.begin());
    std::move(tempStrides.begin(), tempStrides.end(), std::back_inserter(strides));
  }

  // 4. Send the data from all ranks to the master rank
  const auto& boostCommunicator = communicator.getBoostCommunicator();
  const auto numberOfRanks= communicator.size() ;
  std::vector<boost::mpi::request> requests(static_cast<size_t>(numberOfRanks)-1);

  for (int rank = 1; rank < numberOfRanks; ++rank) {
    // Determine where to insert the array
    auto start = totalEvents.data() + strides[rank]-1;
    auto length = numberOfEventsPerRank[rank];
    requests.emplace_back(boostCommunicator.irecv(rank, 1, start, static_cast<int>(length)));
  }
  boost::mpi::wait_all(requests.begin(), requests.end());

  return totalEvents;
}

//ConvertToDistributedMD::BoxStructure* ConvertToDistributedMD::generatePreliminaryBoxStructure(const MDEventList& allEvents) {
//  // To build the box structure we need to:
//  // 1. Create a box controller
//  // 2. Create a seed box
//  // 3. Add the events one by one
//
//  // Create box controller
//
//  //data = new MDBox<MDE, nd>(m_BoxController.get(), 0);
//
//}




  } // namespace MDAlgorithms
} // namespace Mantid
