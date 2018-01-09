#ifndef MANTID_MDALGORITHMS_CONVERTTODISTRIBUTEDMD_H_
#define MANTID_MDALGORITHMS_CONVERTTODISTRIBUTEDMD_H_

#include "MantidMDAlgorithms/DllConfig.h"
#include "MantidMDAlgorithms/BoxStructureSerializer.h"
#include "MantidMDAlgorithms/DistributedCommon.h"
#include "MantidAPI/BoxController.h"
#include "MantidAPI/ParallelAlgorithm.h"
#include "MantidDataObjects/MDEvent.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidDataObjects/MDBoxBase.h"
#include "MantidDataObjects/MDLeanEvent.h"
#include "MantidGeometry/MDGeometry/MDFrameFactory.h"
#include "MantidDataObjects/MDEventFactory.h"
#include "MantidGeometry/Instrument.h"

namespace Mantid {
namespace MDAlgorithms {

struct FrameInformation {
  std::string dimensionNames[3];
  Mantid::Kernel::SpecialCoordinateSystem specialCoordinateSystem;
  Mantid::Geometry::MDFrame_uptr frame;
};


class MANTID_MDALGORITHMS_DLL ConvertToDistributedMD
    : public API::ParallelAlgorithm {
public:
  const std::string name() const override;
  int version() const override;
  const std::string category() const override;
  const std::string summary() const override;

  static const std::string boxSplittingGroupName;

private:
  void init() override;
  void exec() override;

  void setupPreliminaryBoxStructure(DistributedCommon::MDEventList &mdEvents);

  /// Send MDEvents to master
  void sendMDEventsToMaster(const Mantid::Parallel::Communicator &communicator,
                            DistributedCommon::MDEventList &mdEvents) const;

  /// Receive MDEvents on master
  DistributedCommon::MDEventList
  receiveMDEventsOnMaster(const Mantid::Parallel::Communicator &communicator,
                          const DistributedCommon::MDEventList &mdEvents) const;

  DistributedCommon::BoxStructureInformation
  generatePreliminaryBoxStructure(const DistributedCommon::MDEventList &allEvents) const;

  void initBoxControllerProps(const std::string &SplitInto = "5",
                              int SplitThreshold = 1000,
                              int MaxRecursionDepth = 5);

  DistributedCommon::MDEventList
  getFractionEvents(const Mantid::DataObjects::EventWorkspace &workspace, double fraction) const;

  void addEventsToPreliminaryBoxStructure(DistributedCommon::MDEventList& allEvents);

  void redistributeData();

  std::vector<size_t> getBoxPerDepthInformation(const Mantid::Parallel::Communicator &communicator,
                                                std::vector<size_t> numberMDBoxesPerDepth) const;

  std::unordered_map<size_t, std::vector<int>> getRelevantEventsPerRankPerBox(const Mantid::Parallel::Communicator& communicator,
                                                                                   const std::vector<Mantid::API::IMDNode*>& boxes);

  std::unordered_map<size_t, DistributedCommon::MDEventList> sendDataToCorrectRank(const Mantid::Parallel::Communicator& communicator,
                                                                const std::unordered_map<size_t, std::vector<uint64_t>>& relevantEventsPerRankPerBox,
                                                                const std::vector<Mantid::DataObjects::MDBox<Mantid::DataObjects::MDLeanEvent<DistributedCommon::DIM_DISTRIBUTED_TEST>, DistributedCommon::DIM_DISTRIBUTED_TEST>*>& mdBoxes);

  void setNullMDBox(std::vector<Mantid::DataObjects::MDBox<Mantid::DataObjects::MDLeanEvent<DistributedCommon::DIM_DISTRIBUTED_TEST>, DistributedCommon::DIM_DISTRIBUTED_TEST>*>& mdBoxes,
                    size_t startIndex, size_t stopIndex, size_t numberOfDimensions);

  void setRankResponsibilityMap();

  int getResponsibleRank(size_t leafNodeIndex);


    std::vector<coord_t> getWorkspaceExtents() const;

  DistributedCommon::BoxStructureInformation extractBoxStructure(
      Mantid::DataObjects::MDEventWorkspace3Lean &workspace) const;

  void sendRankResponsibility(
      const Mantid::Parallel::Communicator &communicator,
      std::vector<std::pair<size_t, size_t>> &responsibility) const;

  std::vector<std::pair<size_t, size_t>> receiveRankResponsibility(
      const Mantid::Parallel::Communicator &communicator) const;

  SerialBoxStructureInformation serializeBoxStructure(
      const DistributedCommon::BoxStructureInformation &boxStructureInformation) const;

  void sendSerializedBoxStructureInformation(
      const Mantid::Parallel::Communicator &communicator,
      SerialBoxStructureInformation &serializedBoxStructureInformation) const;

  SerialBoxStructureInformation receiveSerializedBoxStructureInformation(
      const Mantid::Parallel::Communicator &communicator) const;

  DistributedCommon::BoxStructureInformation deserializeBoxStructure(
      SerialBoxStructureInformation &serializedBoxStructureInformation) const;

  void continueSplitting();

  void updateMetaData();

  size_t getOffset(int rank, size_t initialMaxID, const std::vector<size_t>& maxIDs) const;

  // -------------------------------------------------------------------------------------------------------------------
  // Methods to build up a temporary MDEventWorkspace
  // -------------------------------------------------------------------------------------------------------------------
  FrameInformation createFrame() const;

  Mantid::DataObjects::MDEventWorkspace3Lean::sptr
  createTemporaryWorkspace() const;

  void addMDEventsToMDEventWorkspace(
      Mantid::DataObjects::MDEventWorkspace3Lean &workspace,
      const DistributedCommon::MDEventList &allEvents) const;

  void setBoxController(Mantid::API::BoxController_sptr bc) const;

  // --------------------------
  // Members
  // --------------------------
  /// The box structure information which contains the box controller and the actuall box structure (of the local rank)
  DistributedCommon::BoxStructureInformation m_boxStructureInformation;
  /// Responsibility map which maps a rank to a start and stop index pair of the box array which is obtained from
  /// IMDNode::getBoxes. This shows which boxes are associated with which rank.
  std::vector<std::pair<size_t, size_t>> m_responsibility;
  /// Is a variation of the m_responsibility vector, which stores the end index of the box range that a rank is
  /// responsible for against the rank which carries that responsibility
  std::unordered_map<size_t, int> m_endBoxIndexRangeVsRank;
  std::vector<size_t> m_endBoxIndexRange;

  size_t m_maxIDBeforeSplit;
};

} // namespace MDAlgorithms
} // namespace Mantid

#endif /* MANTID_MDALGORITHMS_CONVERTTODISTRIBUTEDMD_H_ */
