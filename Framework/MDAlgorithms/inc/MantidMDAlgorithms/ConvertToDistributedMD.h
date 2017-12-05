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
                            const DistributedCommon::MDEventList &mdEvents) const;

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
  DistributedCommon::BoxStructureInformation m_boxStructureInformation;
  std::vector<std::pair<size_t, size_t>> m_responsibility;
};

} // namespace MDAlgorithms
} // namespace Mantid

#endif /* MANTID_MDALGORITHMS_CONVERTTODISTRIBUTEDMD_H_ */
