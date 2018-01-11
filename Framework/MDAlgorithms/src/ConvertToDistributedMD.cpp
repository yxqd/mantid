#include "MantidAPI/IEventWorkspace.h"
#include "MantidAPI/BoxControllerAccess.h"
#include "MantidMDAlgorithms/ConvertToDistributedMD.h"
#include "MantidMDAlgorithms/EventToMDEventConverter.h"
#include "MantidMDAlgorithms/RankResponsibility.h"
#include "MantidParallel/Communicator.h"
#include "MantidParallel/Collectives.h"
#include "MantidDataObjects/NullMDBox.h"
#include "MantidKernel/UnitLabelTypes.h"
#include "MantidKernel/Strings.h"
#include "MantidKernel/ArrayProperty.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/StringTokenizer.h"
#include "MantidKernel/ListValidator.h"

#include <boost/serialization/utility.hpp>
#include <fstream>


namespace {
  class TimerParallel {
  public:
    TimerParallel(const Mantid::Parallel::Communicator& comm) : comm(comm) {
      m_start_cpu_total = std::clock();
      m_start_wall_total = std::chrono::high_resolution_clock::now();
    }

    void start () {
      m_start_cpu = std::clock();
      m_start_wall = std::chrono::high_resolution_clock::now();
    }

    void stop () {
      auto stop_cpu = std::clock();
      auto stop_wall = std::chrono::high_resolution_clock::now();
      m_times_cpu.emplace_back((stop_cpu - m_start_cpu)/CLOCKS_PER_SEC);
      m_times_wall.emplace_back(std::chrono::duration<double>(stop_wall - m_start_wall).count());
    };

    void recordNumEvents(size_t numEvents) {
      m_numEvents = numEvents;
    }

    void dump() {
      char cwd[1024];
      if (getcwd(cwd, sizeof(cwd)) != nullptr) {
        std::string base(cwd);
        if (base.find("scarf") != std::string::npos) {
          fileNameBase = "/home/isisg/scarf672/Mantid2/mpi_test/results/result_";
        } else {
          fileNameBase = "/home/anton/builds";
        }
      }

      // Measure the total time
      auto stop_cpu = std::clock();
      auto stop_wall = std::chrono::high_resolution_clock::now();
      m_times_cpu.emplace_back((stop_cpu - m_start_cpu_total)/CLOCKS_PER_SEC);
      m_times_wall.emplace_back(std::chrono::duration<double>(stop_wall - m_start_wall_total).count());

      // -----------------
      // Receive times
      std::vector<Mantid::Parallel::Request> requests;
      std::vector<double> times_cpu;
      std::vector<double> times_wall;


      if (comm.rank() == 0) {
        const auto numRanks = comm.size();
        const auto size = m_times_cpu.size();
        times_cpu.resize(numRanks*size);
        times_wall.resize(numRanks*size);

        std::copy(m_times_cpu.begin(), m_times_cpu.end(), times_cpu.begin());
        std::copy(m_times_wall.begin(), m_times_wall.end(), times_wall.begin());

        const auto offset = m_times_cpu.size();
        for (int rank = 1; rank < numRanks; ++rank) {
          requests.emplace_back(comm.irecv(rank, 1, times_cpu.data()+offset*rank, static_cast<int>(size)));
          requests.emplace_back(comm.irecv(rank, 2, times_wall.data()+offset*rank, static_cast<int>(size)));
        }
      } else {
        requests.emplace_back(comm.isend(0, 1, m_times_cpu.data(), static_cast<int>(m_times_cpu.size())));
        requests.emplace_back(comm.isend(0, 2, m_times_wall.data(), static_cast<int>(m_times_wall.size())));
      }
      wait_all(requests.begin(), requests.end());
      std::vector<Mantid::Parallel::Request>().swap(requests);

      // -----------------
      // Receive other
      std::vector<size_t> numEvents(static_cast<size_t>(comm.size()));
      if (comm.rank() == 0) {
        Mantid::Parallel::gather(comm, m_numEvents, numEvents, 0);
      } else {
        Mantid::Parallel::gather(comm, m_numEvents, 0);
      }

      // -----------------
      // Save to file
      if (comm.rank() == 0) {
        // Save times
        const auto size = comm.size();
        std::string fileName = fileNameBase + std::to_string(size) + std::string(".txt");
        std::fstream stream;
        stream.open(fileName, std::ios::out | std::ios::app);
        for (auto index=0ul; index < times_cpu.size(); ++index) {
           stream << times_cpu[index] <<","<< times_wall[index] <<"\n";
        }

        // Save other
        saveOther(stream, numEvents);

        stream.close();
      }
    }

  private:
    void saveOther(std::fstream& stream, const std::vector<size_t>& other) {
      for (auto index = 0ul; index < other.size()-1; ++index) {
        stream << other[index] << ",";
      }
      stream << other[other.size()-1] <<"\n";
    }

    std::clock_t m_start_cpu_total;
    std::chrono::system_clock::time_point m_start_wall_total;
    std::clock_t m_start_cpu;
    std::chrono::system_clock::time_point m_start_wall;
    std::vector<double> m_times_cpu;
    std::vector<double> m_times_wall;
    const Mantid::Parallel::Communicator& comm;
    std::string fileNameBase;
    size_t m_numEvents;
  };

}


namespace {
 struct Measurement {
   Measurement(int from, int to, int tag, int length) : from(from), to(to), tag(tag), length(length){}
   int from;
   int to;
   int tag;
   int length;
 };


std::vector<Measurement> sendMeasurement;
std::vector<Measurement> recvMeasurement;

void save(const std::vector<Measurement>& measurement, int rank, const std::string& prefix) {
  char cwd[1024];
  std::string fileNameBase;
  if (getcwd(cwd, sizeof(cwd)) != nullptr) {
    std::string base(cwd);
    if (base.find("scarf") != std::string::npos) {
      fileNameBase = "/home/isisg/scarf672/Mantid2/mpi_test/archive/";
    } else {
      fileNameBase = "/home/anton/builds/Mantid_debug_clion/mpi_test/archive";
    }

    fileNameBase += prefix + std::to_string(rank) + ".txt";

    std::fstream stream;
    stream.open(fileNameBase, std::ios::out | std::ios::app);
    for (auto& e : measurement) {
      stream << prefix <<": " << e.from << "->" << e.to << ", tag " << e.tag << " length " << e.length <<"\n";
    }
    stream.close();
  }
}
}


namespace Mantid {
namespace MDAlgorithms {

using Mantid::Kernel::Direction;
using Mantid::API::WorkspaceProperty;
using DistributedCommon::DIM_DISTRIBUTED_TEST;
using DistributedCommon::BoxStructure;
using DistributedCommon::MDEventList;
using DistributedCommon::BoxStructureInformation;
using namespace Mantid::Kernel;
using namespace Mantid::Geometry;
using namespace Mantid::DataObjects;

const std::string ConvertToDistributedMD::boxSplittingGroupName =
    "Box Splitting Settings";

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(ConvertToDistributedMD)


//----------------------------------------------------------------------------------------------

/// Algorithms name for identification. @see Algorithm::name
const std::string ConvertToDistributedMD::name() const {
  return "ConvertToDistributedMD";
}

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
      Kernel::make_unique<WorkspaceProperty<DataObjects::EventWorkspace>>(
          "InputWorkspace", "", Direction::Input),
      "An input workspace.");

  declareProperty(make_unique<Mantid::Kernel::PropertyWithValue<bool>>(
                      "LorentzCorrection", false, Direction::Input),
                  "Correct the weights of events by multiplying by the Lorentz "
                  "formula: sin(theta)^2 / lambda^4");

  declareProperty(
      make_unique<Mantid::Kernel::PropertyWithValue<double>>("Fraction", 0.01f),
      "Fraction of pulse time that should be used to build the "
      "initial data set.\n");

  std::vector<std::string> propOptions{"Q (lab frame)", "Q (sample frame)",
                                       "HKL"};
  declareProperty(
      "OutputDimensions", "Q (lab frame)",
      boost::make_shared<Mantid::Kernel::StringListValidator>(propOptions),
      "What will be the dimensions of the output workspace?\n"
      "  Q (lab frame): Wave-vector change of the lattice in the lab frame.\n"
      "  Q (sample frame): Wave-vector change of the lattice in the frame of "
      "the sample (taking out goniometer rotation).\n"
      "  HKL: Use the sample's UB matrix to convert to crystal's HKL indices.");

  // ---------------------------------
  // Box Controller settings
  // ---------------------------------
  this->initBoxControllerProps("2" /*SplitInto*/, 1500 /*SplitThreshold*/,
                               20 /*MaxRecursionDepth*/);

  declareProperty(
      make_unique<PropertyWithValue<int>>("MinRecursionDepth", 0),
      "Optional. If specified, then all the boxes will be split to this "
      "minimum recursion depth. 1 = one level of splitting, etc.\n"
      "Be careful using this since it can quickly create a huge number of "
      "boxes = (SplitInto ^ (MinRercursionDepth x NumDimensions)).\n"
      "But setting this property equal to MaxRecursionDepth property is "
      "necessary if one wants to generate multiple file based workspaces in "
      "order to merge them later\n");
  setPropertyGroup("MinRecursionDepth", boxSplittingGroupName);

  std::vector<double> extents{-50., 50., -50., 50., -50., 50.};
  declareProperty(
      Kernel::make_unique<ArrayProperty<double>>("Extents", extents),
      "A comma separated list of min, max for each dimension,\n"
      "specifying the extents of each dimension. Optional, default "
      "+-50 in each dimension.");
  setPropertyGroup("Extents", boxSplittingGroupName);
}

/**
 * TODO: Find a good way to share with BoxSettingsAlgorithm
 */
void ConvertToDistributedMD::initBoxControllerProps(
    const std::string &SplitInto, int SplitThreshold, int MaxRecursionDepth) {
  auto mustBePositive = boost::make_shared<BoundedValidator<int>>();
  mustBePositive->setLower(0);
  auto mustBeMoreThen1 = boost::make_shared<BoundedValidator<int>>();
  mustBeMoreThen1->setLower(1);

  // Split up comma-separated properties
  typedef Mantid::Kernel::StringTokenizer tokenizer;
  tokenizer values(SplitInto, ",",
                   tokenizer::TOK_IGNORE_EMPTY | tokenizer::TOK_TRIM);
  std::vector<int> valueVec;
  valueVec.reserve(values.count());
  for (const auto &value : values)
    valueVec.push_back(boost::lexical_cast<int>(value));

  declareProperty(
      Kernel::make_unique<ArrayProperty<int>>("SplitInto", valueVec),
      "A comma separated list of into how many sub-grid elements each "
      "dimension should split; "
      "or just one to split into the same number for all dimensions. Default " +
          SplitInto + ".");

  declareProperty(
      make_unique<PropertyWithValue<int>>("SplitThreshold", SplitThreshold,
                                          mustBePositive),
      "How many events in a box before it should be split. Default " +
          Kernel::Strings::toString(SplitThreshold) + ".");

  declareProperty(make_unique<PropertyWithValue<int>>(
                      "MaxRecursionDepth", MaxRecursionDepth, mustBeMoreThen1),
                  "How many levels of box splitting recursion are allowed. "
                  "The smallest box will have each side length :math:`l = "
                  "(extents) / (SplitInto^{MaxRecursionDepth}).` "
                  "Default " +
                      Kernel::Strings::toString(MaxRecursionDepth) + ".");
  setPropertyGroup("SplitInto", boxSplittingGroupName);
  setPropertyGroup("SplitThreshold", boxSplittingGroupName);
  setPropertyGroup("MaxRecursionDepth", boxSplittingGroupName);
}

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void ConvertToDistributedMD::exec() {
  // The generation of the distributed MD data requires the following steps:
  // 1. Convert the data which corresponds to the first n percent of the total
  // measurement
  // 2. Build a preliminary box structure on the master rank
  // 3. Determine how the data will be split based on the preliminary box
  // structure and share this information with all
  //    ranks
  // 4. Share the preliminary box structure with all ranks
  // 5. Convert all events
  // 6. Disable the box controller and add events to the local box structure
  // 7. Send data from the each rank to the correct rank
  // 8. Enable the box controller and start splitting the data
  // 9. Ensure that the fileIDs and the box controller stats are correct.
  // 10. Maybe save in this algorithm already

  TimerParallel timer(this->communicator());

  const auto localRank = this->communicator().rank();


  std::cout << "Starting rank " << localRank <<"\n";
  // Get the users's inputs
  EventWorkspace_sptr inputWorkspace = getProperty("InputWorkspace");
  {
    // ----------------------------------------------------------
    // 1. Get a n-percent fraction
    // ----------------------------------------------------------
    double fraction = getProperty("Fraction");
    timer.start();
    auto nPercentEvents = getFractionEvents(*inputWorkspace, fraction);
    timer.stop();

    // -----------------------------------------------------------------
    // 2. + 3. = 4.  Get the preliminary box structure and the partition
    // behaviour
    // -----------------------------------------------------------------
    timer.start();
    setupPreliminaryBoxStructure(nPercentEvents);
    timer.stop();
  }


  // ----------------------------------------------------------
  // 5. Convert all events
  // ----------------------------------------------------------
  {
    timer.start();
    auto allEvents = getFractionEvents(*inputWorkspace, 1.);
    timer.stop();

    // ----------------------------------------------------------
    // 6. Add the local data to the preliminary box structure
    // ----------------------------------------------------------
    timer.start();
    addEventsToPreliminaryBoxStructure(allEvents);
    timer.stop();
  }

 // std::cout << "Finished setting up the box strucutre on rank " << localRank <<"\n";


  // ------------------------------------------------------
  // 7. Redistribute data
  // ----------------------------------------------------------
  timer.start();
  redistributeData();
  timer.stop();
  //std::cout << "Finished redistributing the data on " << localRank <<"\n";
#if 0
  // ----------------------------------------------------------
  // 8. Continue to split locally
  // ----------------------------------------------------------
  timer.start();
  continueSplitting();
  timer.stop();

  // --------------------------------------- -------------------
  // 9. Ensure that box controller and fileIDs are correct
  // ----------------------------------------------------------
  timer.start();
  updateMetaData();
  timer.stop();

  // ----------------------------------------------------------
  // 9. Save?
  // ----------------------------------------------------------
  const auto numEvents = m_boxStructureInformation.boxStructure->getNPoints();
  timer.recordNumEvents(numEvents);
  timer.dump();
#endif
}


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// Definitions for: 1. Convert first n percent of the data
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
MDEventList
ConvertToDistributedMD::getFractionEvents(const EventWorkspace &workspace,
                                          double fraction) const {
  // Get user inputs
  auto extents = getWorkspaceExtents();

  EventToMDEventConverter converter;
  return converter.getEvents(workspace, extents, fraction, QFrame::QLab);
}

/**
 * Fetches the user-specified workspace extents
 * @return a vector with extents (6 values)
 */
std::vector<Mantid::coord_t>
ConvertToDistributedMD::getWorkspaceExtents() const {
  std::vector<double> extentsInput = getProperty("Extents");
  // Replicate a single min,max into several
  if (extentsInput.size() != DIM_DISTRIBUTED_TEST * 2)
    throw std::invalid_argument(
      "You must specify multiple of 2 extents, ie two for each dimension.");

  // Convert input to coordinate type
  std::vector<coord_t> extents;
  std::transform(extentsInput.begin(), extentsInput.end(),
                 std::back_inserter(extents),
                 [](double value) { return static_cast<coord_t>(value); });

  return extents;
}

// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// Definitions for: 2 + 3 + 4. Get the prelim. box structure
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
void ConvertToDistributedMD::setupPreliminaryBoxStructure(
  MDEventList &mdEvents) {
  // The following steps are performed
  // 1. The first data is sent from all the ranks to the master rank
  // 2. The box structure is built on the master rank
  // 3. The responsibility of the ranks is determined, it is decided which ranks
  // will take care of which boxes
  // 4. This responsibility is shared with all ranks.
  // 5. The master rank serializes the box structure
  // 6. The serialized box structure is broadcast from the master rank to all
  // other ranks
  // 7. The other ranks deserialize the box structure (but without the event
  // information)
  const auto &communicator = this->communicator();

  SerialBoxStructureInformation serialBoxStructureInformation;

  if (communicator.rank() == 0) {
    // 1.b Receive the data from all the other ranks
    auto allEvents = receiveMDEventsOnMaster(communicator, mdEvents);

    // 2. Build the box structure on the master rank
    auto boxStructureInformation = generatePreliminaryBoxStructure(allEvents);

    // 3. Determine splitting
    RankResponsibility rankResponsibility;
    m_responsibility = rankResponsibility.getResponsibilites(
      communicator.size(), boxStructureInformation.boxStructure.get());

    // 4.a Broadcast splitting behaviour
    sendRankResponsibility(communicator, m_responsibility);

    // 5. Serialize the box structure
    serialBoxStructureInformation =
      serializeBoxStructure(boxStructureInformation);

    // 6.a Broadcast serialized box structure to all other ranks
    sendSerializedBoxStructureInformation(communicator,
                                          serialBoxStructureInformation);

  } else {
    // 1.a Send the event data to the master rank
    sendMDEventsToMaster(communicator, mdEvents);

    // 4.b Receive the splitting behaviour from the broadcast
    m_responsibility = receiveRankResponsibility(communicator);

    // 6.b Receive the box structure from master
    serialBoxStructureInformation =
      receiveSerializedBoxStructureInformation(communicator);
  }

  // 7. Deserialize on all ranks
  auto hollowBoxStructureInformation =
    deserializeBoxStructure(serialBoxStructureInformation);

  // Create a rank responsibility map
  setRankResponsibilityMap();


  // Set results
  m_boxStructureInformation = std::move(hollowBoxStructureInformation);
}

void ConvertToDistributedMD::setRankResponsibilityMap() {
  for (auto rank = 0ul; rank < m_responsibility.size(); ++rank) {
    m_endBoxIndexRangeVsRank.emplace(m_responsibility[rank].second, rank);
    m_endBoxIndexRange.push_back(m_responsibility[rank].second);
  }
}


void ConvertToDistributedMD::sendMDEventsToMaster(
  const Mantid::Parallel::Communicator &communicator,
  MDEventList &mdEvents) const {
  // Send the totalNumberOfEvents
  auto totalNumberEvents = mdEvents.size();
  gather(communicator, totalNumberEvents, 0);

  // Send the vector of md events
  auto sizeMdEvents = sizeof(std::remove_reference<decltype(mdEvents)>::type::value_type);
  communicator.send(0, 2, reinterpret_cast<char*>(mdEvents.data()),
                    static_cast<int>(mdEvents.size()*sizeMdEvents));
}


MDEventList ConvertToDistributedMD::receiveMDEventsOnMaster(
  const Mantid::Parallel::Communicator &communicator,
  const MDEventList &mdEvents) const {
  // In order to get all relevant md events onto the master rank we need to:
  // 1. Inform the master rank about the number of events on all ranks (gather)
  // 2. Create a sufficiently large buff on master
  // 3. Determine the stride that each rank requires, ie which offset each rank
  // requires
  // 4. Send the data to master. This would be a natural operation for gatherv,
  // however this is not available in
  //    boost 1.58 (in 1.59+ it is).

  // 1. Get the totalNumberOfEvents for each rank
  std::vector<size_t> numberOfEventsPerRank;
  auto masterNumberOfEvents = mdEvents.size();
  gather(communicator, masterNumberOfEvents, numberOfEventsPerRank, 0);

  // 2. Create a buffer
  auto totalNumberOfEvents = std::accumulate(numberOfEventsPerRank.begin(),
                                             numberOfEventsPerRank.end(), 0ul);
  MDEventList totalEvents(totalNumberOfEvents);
  std::copy(mdEvents.begin(), mdEvents.end(), totalEvents.begin());

  // 3. Determine the strides
  std::vector<int> strides;
  strides.reserve(numberOfEventsPerRank.size());
  strides.push_back(0);
  {
    std::vector<int> tempStrides(numberOfEventsPerRank.size() - 1);
    std::partial_sum(numberOfEventsPerRank.begin(),
                     numberOfEventsPerRank.end() - 1, tempStrides.begin());
    std::move(tempStrides.begin(), tempStrides.end(),
              std::back_inserter(strides));
  }

  // 4. Send the data from all ranks to the master rank
  const auto numberOfRanks = communicator.size();
  std::vector<Mantid::Parallel::Request> requests;
  requests.reserve(static_cast<size_t>(numberOfRanks) - 1);
  auto sizeMdEvents = sizeof(std::remove_reference<decltype(mdEvents)>::type::value_type);
  for (int rank = 1; rank < numberOfRanks; ++rank) {
    // Determine where to insert the array
    auto start = totalEvents.data() + strides[rank];
    auto length = numberOfEventsPerRank[rank];
    requests.emplace_back(
      communicator.irecv(rank, 2, reinterpret_cast<char*>(start), static_cast<int>(length*sizeMdEvents)));
  }
  wait_all(requests.begin(), requests.end());
  return totalEvents;
}


BoxStructureInformation ConvertToDistributedMD::generatePreliminaryBoxStructure(
  const MDEventList &allEvents) const {
  // To build the box structure we need to:
  // 1. Generate a temporary MDEventWorkspace locally on the master rank
  // 2. Populate the workspace with the MDEvents
  // 3. Extract the box structure and the box controller from the workspace, ie
  //    remove ownership from the underlying data structures

  // 1. Create temporary MDEventWorkspace
  auto tempWorkspace = createTemporaryWorkspace();

  // 2. Populate workspace with MDEvents
  addMDEventsToMDEventWorkspace(*tempWorkspace, allEvents);

  // 3. Extract box structure and box controller
  return extractBoxStructure(*tempWorkspace);
}

void ConvertToDistributedMD::sendSerializedBoxStructureInformation(
  const Mantid::Parallel::Communicator &communicator,
  SerialBoxStructureInformation &serializedBoxStructureInformation) const {
  boost::mpi::broadcast(communicator, serializedBoxStructureInformation, 0);
}

SerialBoxStructureInformation
ConvertToDistributedMD::receiveSerializedBoxStructureInformation(
  const Mantid::Parallel::Communicator &communicator) const {
  SerialBoxStructureInformation serializedBoxStructureInformation;
  boost::mpi::broadcast(communicator, serializedBoxStructureInformation, 0);
  return serializedBoxStructureInformation;
}

void ConvertToDistributedMD::sendRankResponsibility(
  const Mantid::Parallel::Communicator &communicator,
  std::vector<std::pair<size_t, size_t>> &responsibility) const {
  boost::mpi::broadcast(communicator, responsibility.data(),
                        communicator.size(), 0);
}

std::vector<std::pair<size_t, size_t>>
ConvertToDistributedMD::receiveRankResponsibility(
  const Mantid::Parallel::Communicator &communicator) const {
  std::vector<std::pair<size_t, size_t>> responsibility;
  responsibility.resize(static_cast<size_t>(communicator.size()));
  boost::mpi::broadcast(communicator, responsibility.data(),
                        communicator.size(), 0);
  return responsibility;
}

SerialBoxStructureInformation ConvertToDistributedMD::serializeBoxStructure(
  const BoxStructureInformation &boxStructureInformation) const {
  BoxStructureSerializer serializer;
  return serializer.serializeBoxStructure(boxStructureInformation);
}


BoxStructureInformation ConvertToDistributedMD::deserializeBoxStructure(
  SerialBoxStructureInformation &serializedBoxStructureInformation) const {
  BoxStructureSerializer serializer;
  return serializer.deserializeBoxStructure(serializedBoxStructureInformation);
}


void ConvertToDistributedMD::addMDEventsToMDEventWorkspace(
  MDEventWorkspace3Lean &workspace, const MDEventList &allEvents) const {
  // We could call the addEvents method, but it seems to behave slightly
  // differently
  // to addEvent. TODO: Investigate if addEvents could be used.
  size_t eventsAddedSinceLastSplit = 0;
  size_t eventsAddedTotal = 0;
  auto boxController = workspace.getBoxController();
  size_t lastNumBoxes = boxController->getTotalNumMDBoxes();

  for (const auto &event : allEvents) {
    // We have a different situation here compared to ConvertToDiffractionMD,
    // since we already have
    // all events converted
    if (boxController->shouldSplitBoxes(
      eventsAddedTotal, eventsAddedSinceLastSplit, lastNumBoxes)) {
      workspace.splitAllIfNeeded(nullptr);
      eventsAddedSinceLastSplit = 0;
      lastNumBoxes = boxController->getTotalNumMDBoxes();
    }
    ++eventsAddedSinceLastSplit;
    ++eventsAddedTotal;
    workspace.addEvent(event);
  }
  workspace.splitAllIfNeeded(nullptr);
  workspace.refreshCache();
}

BoxStructureInformation ConvertToDistributedMD::extractBoxStructure(
  MDEventWorkspace3Lean &workspace) const {
  // Extract the box controller
  BoxStructureInformation boxStructureInformation;
  workspace.transferInternals(boxStructureInformation.boxController,
                              boxStructureInformation.boxStructure);
  return boxStructureInformation;
}


MDEventWorkspace3Lean::sptr
ConvertToDistributedMD::createTemporaryWorkspace() const {
  auto imdWorkspace = DataObjects::MDEventFactory::CreateMDWorkspace(
    DIM_DISTRIBUTED_TEST, "MDLeanEvent");
  auto workspace =
    boost::dynamic_pointer_cast<DataObjects::MDEventWorkspace3Lean>(
      imdWorkspace);
  auto frameInformation = createFrame();
  auto extents = getWorkspaceExtents();
  for (size_t d = 0; d < DIM_DISTRIBUTED_TEST; d++) {
    MDHistoDimension *dim = new MDHistoDimension(
      frameInformation.dimensionNames[d], frameInformation.dimensionNames[d],
      *frameInformation.frame, static_cast<coord_t>(extents[d * 2]),
      static_cast<coord_t>(extents[d * 2 + 1]), 10);
    workspace->addDimension(MDHistoDimension_sptr(dim));
  }
  workspace->initialize();

  // Build up the box controller
  auto bc = workspace->getBoxController();
  this->setBoxController(bc);
  workspace->splitBox();

  // Perform minimum recursion depth splitting
  int minDepth = this->getProperty("MinRecursionDepth");
  int maxDepth = this->getProperty("MaxRecursionDepth");
  if (minDepth > maxDepth)
    throw std::invalid_argument(
      "MinRecursionDepth must be <= MaxRecursionDepth ");
  workspace->setMinRecursionDepth(size_t(minDepth));

  workspace->setCoordinateSystem(frameInformation.specialCoordinateSystem);
  return workspace;
}

FrameInformation ConvertToDistributedMD::createFrame() const {
  using Mantid::Kernel::SpecialCoordinateSystem;
  FrameInformation information;
  std::string outputDimensions = getProperty("OutputDimensions");
  auto frameFactory = makeMDFrameFactoryChain();

  if (outputDimensions == "Q (sample frame)") {
    // Names
    information.dimensionNames[0] = "Q_sample_x";
    information.dimensionNames[1] = "Q_sample_y";
    information.dimensionNames[2] = "Q_sample_z";
    information.specialCoordinateSystem = Mantid::Kernel::QSample;
    // Frame
    MDFrameArgument frameArgQSample(QSample::QSampleName, "");
    information.frame = frameFactory->create(frameArgQSample);
  } else if (outputDimensions == "HKL") {
    information.dimensionNames[0] = "H";
    information.dimensionNames[1] = "K";
    information.dimensionNames[2] = "L";
    information.specialCoordinateSystem = Mantid::Kernel::HKL;
    MDFrameArgument frameArgQLab(HKL::HKLName,
                                 Mantid::Kernel::Units::Symbol::RLU.ascii());
    information.frame = frameFactory->create(frameArgQLab);
  } else {
    information.dimensionNames[0] = "Q_lab_x";
    information.dimensionNames[1] = "Q_lab_y";
    information.dimensionNames[2] = "Q_lab_z";
    information.specialCoordinateSystem = Mantid::Kernel::QLab;
    MDFrameArgument frameArgQLab(QLab::QLabName, "");
    information.frame = frameFactory->create(frameArgQLab);
  }

  return information;
}

void ConvertToDistributedMD::setBoxController(
  Mantid::API::BoxController_sptr bc) const {
  size_t numberOfDimensions = bc->getNDims();

  int splitThreshold = this->getProperty("SplitThreshold");
  bc->setSplitThreshold(static_cast<size_t>(splitThreshold));
  int maxRecursionDepth = this->getProperty("MaxRecursionDepth");
  bc->setMaxDepth(static_cast<size_t>(maxRecursionDepth));

  std::vector<int> splits = getProperty("SplitInto");
  if (splits.size() == 1) {
    bc->setSplitInto(static_cast<size_t>(splits[0]));
  } else if (splits.size() == numberOfDimensions) {
    for (size_t d = 0; d < numberOfDimensions; ++d)
      bc->setSplitInto(d, static_cast<size_t>(splits[d]));
  } else
    throw std::invalid_argument("SplitInto parameter has " +
                                Strings::toString(splits.size()) +
                                " arguments. It should have either 1, or the "
                                  "same as the number of dimensions.");
  bc->resetNumBoxes();
}


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// Definitions for: 5 + 6 Add all events to the box structure
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
void ConvertToDistributedMD::addEventsToPreliminaryBoxStructure(
  DistributedCommon::MDEventList &allEvents) {
  // Add all events to the hollow box structure
  for (auto &event : allEvents) {
    // TODO: create a move-enabled addEvent
    m_boxStructureInformation.boxStructure->addEvent(std::move(event));
  }
}


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// Definitions for: 7. Redistribute data
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

void ConvertToDistributedMD::redistributeData() {
  const auto &communicator = this->communicator();
  // Build up the data which needs to be transmitted
  std::vector<Mantid::API::IMDNode *> boxes;
  m_boxStructureInformation.boxStructure->getBoxes(boxes, 1000, true);
  std::vector<MDBox<MDLeanEvent<DIM_DISTRIBUTED_TEST>, DIM_DISTRIBUTED_TEST> *>
    mdBoxes;
  for (auto &box : boxes) {
    mdBoxes.emplace_back(dynamic_cast<
                           MDBox<MDLeanEvent<DIM_DISTRIBUTED_TEST>, DIM_DISTRIBUTED_TEST> *>(box));
  }


  // Determine the number of events per rank per box
  auto relevantEventsPerRankPerBox =
    getRelevantEventsPerRankPerBox(communicator, boxes);

  // Send the actual data
  auto boxVsMDEvents =
    sendDataToCorrectRank(communicator, relevantEventsPerRankPerBox, mdBoxes);

  // Place the data into the correct boxes
  auto localRank = communicator.rank();
  auto startIndex = m_responsibility[localRank].first;
  auto stopIndex = m_responsibility[localRank].second;
  for (auto index = startIndex; index <= stopIndex; ++index) {
    auto mdBox = mdBoxes[index];

    auto &events = boxVsMDEvents.at(index);
    auto& oldEvents = mdBox->getEvents();
    oldEvents.swap(events);
    MDEventList().swap(events);
  }

  // Remove the data which is not relevant for the local rank which is before
  // the startIndex
  auto const numberOfDimensions = mdBoxes[0]->getNumDims();
  setNullMDBox(mdBoxes, 0ul, startIndex, numberOfDimensions);

  // Remove the data which is not relevant for the local rank which is after the
  // stopIndex
  setNullMDBox(mdBoxes, stopIndex + 1, boxes.size(), numberOfDimensions);

  // We need to refresh the cache
  m_boxStructureInformation.boxStructure->refreshCache(nullptr);


  // Cache the max ID, we need it later when updating the fileIDs on all the
  // ranks. Note that getMaxId will return the next available id, ie it is
  // not really the max id.
  m_maxIDBeforeSplit = m_boxStructureInformation.boxController->getMaxId() - 1;
}


void ConvertToDistributedMD::setNullMDBox(std::vector<MDBox<MDLeanEvent<DIM_DISTRIBUTED_TEST>, DIM_DISTRIBUTED_TEST>*>& mdBoxes,
                  size_t startIndex, size_t stopIndex, const size_t numberOfDimensions) {
  for (auto index = startIndex; index < stopIndex; ++index) {
    auto mdBox = mdBoxes[index];

    // On this rank the box has to be an MDEventBox
    if (!mdBox->isBox()) {
      throw std::runtime_error("Expected an MDBox, but got an MDGridBox");
    }

    const auto& boxController  = mdBox->getBoxController();
    const auto depth = mdBox->getDepth();

    std::vector<Mantid::Geometry::MDDimensionExtents<coord_t>> extents;
    for (auto dimension = 0ul; dimension < numberOfDimensions; ++dimension) {
      extents.emplace_back(mdBox->getExtents(dimension));
    }

    const auto boxID = mdBox->getID();
    auto responsibleRank = getResponsibleRank(index);

    // Create a new NullMDBox
    auto parent = mdBox->getParent();
    auto newNullMDBox = Mantid::Kernel::make_unique<NullMDBox<MDLeanEvent<DIM_DISTRIBUTED_TEST>, DIM_DISTRIBUTED_TEST>>(boxController, depth, extents, boxID, responsibleRank);
    newNullMDBox->setShowGlobalValues(false);
    newNullMDBox->setParent(parent);

    // Find the child index and replace it in the parent
    auto mdGridBoxParent = dynamic_cast<MDGridBox<MDLeanEvent<DIM_DISTRIBUTED_TEST>, DIM_DISTRIBUTED_TEST> *>(parent);
    const auto childIndex = mdGridBoxParent->getChildIndexFromID(boxID);
    mdGridBoxParent->setChild(childIndex, newNullMDBox.release());
  }
}


int ConvertToDistributedMD::getResponsibleRank(size_t leafNodeIndex) {
  auto result = std::lower_bound(m_endBoxIndexRange.begin(), m_endBoxIndexRange.end(), leafNodeIndex);
  return m_endBoxIndexRangeVsRank[*result];
}


std::unordered_map<size_t, std::vector<uint64_t>>
ConvertToDistributedMD::getRelevantEventsPerRankPerBox(
  const Mantid::Parallel::Communicator &communicator,
  const std::vector<Mantid::API::IMDNode *> &boxes) {
  // We want to get the number of events per rank per box for the boxes which
  // are relevant for the local
  // rank, ie the boxes for which the local rank is responsible.

  const auto localRank = communicator.rank();
  const auto numberOfRanks = static_cast<size_t>(communicator.size());

  // Initialize Events
  std::unordered_map<size_t, std::vector<uint64_t>> relevantEventsPerRankPerBox;
  auto range = m_responsibility[localRank];
  for (size_t index = range.first; index <= range.second; ++index) {
    auto &elements = relevantEventsPerRankPerBox[index];
    elements.resize(numberOfRanks, 0);
  }

  // Get the ranges for rank 0
  auto rankOfCurrentIndex = 0;
  auto startIndex = m_responsibility[rankOfCurrentIndex].first;
  auto stopIndex = m_responsibility[rankOfCurrentIndex].second;
  auto isInBounds = [&startIndex, &stopIndex](size_t index) {
    return (startIndex <= index) && (index <= stopIndex);
  };

  std::vector<Mantid::Parallel::Request> requests;
  std::vector<uint64_t> nEventBuffer;
  const boost::mpi::communicator& boostComm = communicator;
  const MPI_Comm comm = boostComm;

  std::vector<MPI_Request> mpi_requests;
  std::vector<MPI_Status> mpi_status;

  for (auto index = 0ul; index < boxes.size(); ++index) {

    // Check if we are still on the same rank. If not then we need to
    // change the current rank and the corresponding indices
    if (!isInBounds(index)) {
      ++rankOfCurrentIndex;
      startIndex = m_responsibility[rankOfCurrentIndex].first;
      stopIndex = m_responsibility[rankOfCurrentIndex].second;
    }

    auto box = boxes[index];

    // Now check if the local rank sends or receives a value for this
    // We send a value if the local rank is not the same as the rank
    // which is associated with the box, else we receive a value
    // from all other boxes
    if (localRank == rankOfCurrentIndex) {
      for (auto rank = 0; rank < static_cast<int>(numberOfRanks); ++rank) {
        auto &eventsPerBox = relevantEventsPerRankPerBox.at(index);

        // We don't send anything to ourselves, however we want to record the
        // number of events
        // that we have already in the right place
        if (rank == localRank) {
          eventsPerBox[rank] = box->getNPoints();
          continue;
        }
        mpi_status.emplace_back();
        mpi_requests.emplace_back();
        MPI_Irecv(&eventsPerBox[rank], 1, MPI_UINT64_T, rank, static_cast<int>(index), comm, &mpi_requests.back());
      }
    } else {
      nEventBuffer.emplace_back(box->getNPoints());
      mpi_status.emplace_back();
      mpi_requests.emplace_back();
      MPI_Isend(&nEventBuffer.back(), 1, MPI_UINT64_T, rankOfCurrentIndex, static_cast<int>(index), comm, &mpi_requests.back());
    }
  }

  auto result = MPI_Waitall(static_cast<int>(mpi_requests.size()), mpi_requests.data(), mpi_status.data());
  if (result != MPI_SUCCESS) {
    std::string message = "There has been an issue sharing the event numbers on rank " + std::to_string(localRank);
    throw std::runtime_error(message);
  }

  auto sync = MPI_Barrier(comm);
  if (sync != MPI_SUCCESS) {
    throw std::runtime_error("Sync failed");
  }

  return relevantEventsPerRankPerBox;
}


std::unordered_map<size_t, DistributedCommon::MDEventList>
ConvertToDistributedMD::sendDataToCorrectRank(
  const Mantid::Parallel::Communicator &communicator,
  const std::unordered_map<size_t, std::vector<uint64_t>> &
  relevantEventsPerRankPerBox,
  const std::vector<MDBox<MDLeanEvent<DIM_DISTRIBUTED_TEST>,
    DIM_DISTRIBUTED_TEST> *> &mdBoxes) {
  // We send the data between all nodes. For this we:
  // 1. Set up a buffer which will accumulate the data. The buffer is a map from
  // boxIndex to an event vector
  // 2. Determine the strides. Since we accumulate data for one box from several
  // ranks into a single vector
  //    we need to determine at which starting position of the array the data of
  //    a particular rank should be added.
  // 3. Iterate over all boxes and determine if a box is asociated with the
  // local rank, in which case we want to
  //    receive data from all other ranks, else we need to send it

  // --------------------
  // 1. Set up the buffer
  // --------------------
  const auto localRank = communicator.rank();
  const auto numberOfRanks = static_cast<size_t>(communicator.size());

  std::unordered_map<size_t, MDEventList> boxVsMDEvents;
  auto startIndex = m_responsibility[localRank].first;
  auto stopIndex = m_responsibility[localRank].second;
  for (auto index = startIndex; index <= stopIndex; ++index) {
    // Directly with MDLeanEvents
    auto &numberOfEvents = relevantEventsPerRankPerBox.at(index);
    auto totalNumEvents =
      std::accumulate(numberOfEvents.begin(), numberOfEvents.end(), 0ul);
    auto &entry = boxVsMDEvents[index];
    entry.resize(totalNumEvents);
  }

  if (localRank == 0) {

  }


  // -------------------------
  // 2. Determine the strides
  // -------------------------
  std::unordered_map<size_t, std::vector<uint64_t>> strides;
  for (auto index = startIndex; index <= stopIndex; ++index) {
    auto &stride = strides[index];
    auto &numberOfEventsPerRank = relevantEventsPerRankPerBox.at(index);
    stride.push_back(0);
    std::vector<uint64_t> tempStrides(numberOfEventsPerRank.size() - 1);
    std::partial_sum(numberOfEventsPerRank.begin(),
                     numberOfEventsPerRank.end() - 1, tempStrides.begin());
    std::move(tempStrides.begin(), tempStrides.end(),
              std::back_inserter(stride));
  }


  // -------------------------
  // 3. Send the data
  // -------------------------
  auto rankOfCurrentIndex = 0;
  startIndex = m_responsibility[rankOfCurrentIndex].first;
  stopIndex = m_responsibility[rankOfCurrentIndex].second;
  auto isInBounds = [&startIndex, &stopIndex](size_t index) {
    return (startIndex <= index) && (index <= stopIndex);
  };

  std::vector<Mantid::Parallel::Request> requests;

  const boost::mpi::communicator& boostComm = communicator;
  MPI_Comm comm = boostComm;

  std::vector<MPI_Request> mpi_requests;
  std::vector<MPI_Status> mpi_status;

  for (auto index = 0ul; index < mdBoxes.size(); ++index) {

    // Check if we are still on the same rank. If not then we need to
    // change the current rank and the corresponding indices
    if (!isInBounds(index)) {
      ++rankOfCurrentIndex;
      startIndex = m_responsibility[rankOfCurrentIndex].first;
      stopIndex = m_responsibility[rankOfCurrentIndex].second;
    }

    auto mdBox = mdBoxes[index];

    // Now check if the local rank sends or receives a value for this
    // We send a value if the local rank is not the same as the rank
    // which is associated with the box, else we receive a value
    // from all other boxes
    const auto sizeOfMDLeanEvent = sizeof(Mantid::DataObjects::MDLeanEvent<DIM_DISTRIBUTED_TEST>);
    if (localRank == rankOfCurrentIndex) {
      for (auto rank = 0; rank < static_cast<int>(numberOfRanks); ++rank) {

        const auto &relevantEventsPerRank =
          relevantEventsPerRankPerBox.at(index);
        const auto numberOfEvents = relevantEventsPerRank[rank];

        const auto &stride = strides.at(index);
        const auto offset = stride[rank];

        auto &events = boxVsMDEvents.at(index);

        // If we don't have events, then we don't do anything
        if (numberOfEvents == 0) {
          continue;
        }

        auto insertionPoint = events.data() + offset;

        // We don't send anything to ourselves, however we want to record the
        // number of events
        // that we have already in the right place
        if (rank == localRank) {
          auto& boxEvent = mdBox->getEvents();
          auto start = boxEvent.data();
          auto size = mdBox->getDataInMemorySize();

          if (size != numberOfEvents) {
            throw std::runtime_error("Mismatch in the number of events.");
          }

          for (auto eventIndex = 0ul; eventIndex < size; ++eventIndex) {
            *(insertionPoint + eventIndex) = *(start + eventIndex);
          }
          continue;
        }
        #if 0
        requests.emplace_back(
          communicator.irecv(rank, static_cast<int>(index), reinterpret_cast<char*>(insertionPoint),
                             static_cast<int>(numberOfEvents*sizeOfMDLeanEvent)));
        #else
        mpi_requests.emplace_back();
        mpi_status.emplace_back();
        //std::cout << "RECV: " << rank << " -> " << localRank << ", tag " << index <<" length " << static_cast<int>(numberOfEvents*sizeOfMDLeanEvent) <<"\n";
        MPI_Irecv(reinterpret_cast<char*>(insertionPoint),
                  static_cast<int>(numberOfEvents*sizeOfMDLeanEvent),
                  MPI_CHAR,
                  rank,
                  static_cast<int>(index),
                  comm,
                  &mpi_requests.back());
        recvMeasurement.emplace_back(rank, localRank, index, static_cast<int>(numberOfEvents*sizeOfMDLeanEvent));
        #endif
      }
    } else {
      auto& events = mdBox->getEvents();
      if (!events.empty()) {
        #if 0
        requests.emplace_back(communicator.isend(
          rankOfCurrentIndex, static_cast<int>(index), reinterpret_cast<char*>(events.data()),
          static_cast<int>(mdBox->getDataInMemorySize()*sizeOfMDLeanEvent)));
        #else
        mpi_requests.emplace_back();
        mpi_status.emplace_back();
       // std::cout << "SEND: " << localRank << " -> " << rankOfCurrentIndex << ", tag " << index <<" length " << static_cast<int>(mdBox->getDataInMemorySize()*sizeOfMDLeanEvent) <<"\n";
        MPI_Isend(reinterpret_cast<char*>(events.data()),
                  static_cast<int>(mdBox->getDataInMemorySize()*sizeOfMDLeanEvent),
                  MPI_CHAR,
                  rankOfCurrentIndex,
                  static_cast<int>(index),
                  comm,
                  &mpi_requests.back());
        sendMeasurement.emplace_back(localRank, rankOfCurrentIndex, index, static_cast<int>(mdBox->getDataInMemorySize()*sizeOfMDLeanEvent));
      #endif
      }
    }
  }

#if 0
  Mantid::Parallel::wait_all(requests.begin(), requests.end());
#else
  auto result = MPI_Waitall(static_cast<int>(mpi_requests.size()), mpi_requests.data(), mpi_status.data());
  if (result != MPI_SUCCESS) {
    std::string message = "There has been an issue sharing the event numbers on rank " + std::to_string(localRank);
    throw std::runtime_error(message);
  } else {
    //std::cout << "Finished sending data on rank " << localRank <<"\n";
  }
#endif
  auto sync = MPI_Barrier(comm);
  if (sync != MPI_SUCCESS) {
    throw std::runtime_error("Sync failed");
  }

  //std::cout << "GOT HERE " << localRank <<"\n";
  save(sendMeasurement, localRank, "SEND");
  save(recvMeasurement, localRank, "RECV");

  return boxVsMDEvents;
}


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// Definitions for: 8. Continue splitting
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
void ConvertToDistributedMD::continueSplitting() {
  auto root = m_boxStructureInformation.boxStructure.get();
  root->splitAllIfNeeded(nullptr);
}


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// Definitions for: 9. Update meta data
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

void ConvertToDistributedMD::updateMetaData() {
  // Since the box controller was duplicated onto several ranks we have several
  // things that we need to correct
  // 1. The fileID in the boxes will be incorrect for most boxes.
  // 2. A box controller which has the information about all of the boxes needs
  // to be created
  // 3. Update the data in the NullMDBoxes

  // --------------------------
  // 1. Updating file IDs
  //   a. Get the max file ID on the each rank and share
  //   b. Have each rank update its file IDs based on the file ID range of the
  //   other ranks
  // --------------------------
  const auto &communicator = this->communicator();
  std::vector<size_t> maxIds;
  all_gather(communicator, m_boxStructureInformation.boxController->getMaxId()-1,
             maxIds);

  const auto localRank = communicator.rank();
  const auto offset = getOffset(localRank, m_maxIDBeforeSplit, maxIds);
  std::vector<API::IMDNode *> boxes;
  m_boxStructureInformation.boxStructure->getBoxes(boxes, 1000,
                                                   false /*leaves only*/);
  for (auto box : boxes) {
    // When only need to update the boxID if it is part of a newly created box
    auto boxID = box->getID();
    if (boxID > m_maxIDBeforeSplit) {
      boxID += offset;
      box->setID(boxID);
    }
  }

  // --------------------------------------
  // 2. Update the box controller settings. This requires updating:
  //   a. The maxID
  //   b. Number of boxes per level
  //   c. Number of grid boxes per level
  //   d. The maximum number of boxes per level / this can be reset via a method
  // --------------------------------------
  // MaxID
  const auto offsetLastRank =
    getOffset(communicator.size() - 1, m_maxIDBeforeSplit, maxIds);

  const auto maxId = maxIds[communicator.size() - 1] + offsetLastRank;
  // Number of boxes per level
  const auto mdBoxPerDepth = getBoxPerDepthInformation(
    communicator, m_boxStructureInformation.boxController->getNumMDBoxes());

  // Number of grid boxes per level
  const auto mdGridBoxPerDepth = getBoxPerDepthInformation(
    communicator,
    m_boxStructureInformation.boxController->getNumMDGridBoxes());

  // Apply changes to the box controller
  Mantid::API::BoxControllerAccess access;
  auto &boxController = *m_boxStructureInformation.boxController;
  boxController.setMaxId(maxId);
  access.setNumMDBoxes(boxController, mdBoxPerDepth);
  access.setNumMDGridBoxes(boxController, mdGridBoxPerDepth);


  // --------------------------------------
  // 2. Update values in NullMDBoxes
  //    TODO: Update the signals stored in the NullMDBoxes. This would be nice to have as a check but
  //          to check scaling up to this point.
  // --------------------------------------
}


size_t
ConvertToDistributedMD::getOffset(int rank, size_t initialMaxID,
                                  const std::vector<size_t> &maxIDs) const {
  // Calculate what the number of new boxes on each rank is
  auto numberOfBoxesOnRanks = maxIDs;
  std::transform(numberOfBoxesOnRanks.begin(),
                 numberOfBoxesOnRanks.begin() + rank,
                 numberOfBoxesOnRanks.begin(),
                 [initialMaxID](size_t value) { return value - initialMaxID; });
  return std::accumulate(numberOfBoxesOnRanks.begin(),
                         numberOfBoxesOnRanks.begin() + rank, 0ul);
}


std::vector<size_t> ConvertToDistributedMD::getBoxPerDepthInformation(
    const Mantid::Parallel::Communicator &communicator,
    std::vector<size_t> numberBoxesPerDepth) const {
  std::vector<size_t> depthBoxes;
  all_gather(communicator, numberBoxesPerDepth.size(), depthBoxes);
  const auto maxDepth = *std::max_element(depthBoxes.begin(), depthBoxes.end());

  numberBoxesPerDepth.resize(maxDepth);
  std::vector<size_t> allNumberBoxes;
  all_gather(communicator, numberBoxesPerDepth.data(),
             static_cast<int>(numberBoxesPerDepth.size()), allNumberBoxes);

  std::vector<size_t> accumulatedAllNumberBoxes;
  const auto numberOfRanks = static_cast<size_t>(communicator.size());
  for (auto depth = 0ul; depth < maxDepth; ++depth) {
    auto sum = 0ul;
    for (auto rank = 0ul; rank < numberOfRanks; ++rank) {
      sum += allNumberBoxes[depth + rank * maxDepth];
    }
    accumulatedAllNumberBoxes.push_back(sum);
  }
  return accumulatedAllNumberBoxes;
}

} // namespace MDAlgorithms
} // namespace Mantid
