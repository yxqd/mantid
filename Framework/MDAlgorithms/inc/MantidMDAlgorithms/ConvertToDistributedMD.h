#ifndef MANTID_MDALGORITHMS_CONVERTTODISTRIBUTEDMD_H_
#define MANTID_MDALGORITHMS_CONVERTTODISTRIBUTEDMD_H_

#include "MantidMDAlgorithms/DllConfig.h"
#include "MantidAPI/ParallelAlgorithm.h"
#include "MantidDataObjects/MDEvent.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidDataObjects/MDBoxBase.h"
#include "MantidDataObjects/MDLeanEvent.h"
#include "MantidAPI/ParallelAlgorithm.h"
#include "MantidGeometry/MDGeometry/MDFrameFactory.h"
#include "MantidDataObjects/MDEventFactory.h"
#include "MantidGeometry/Instrument.h"
#include "MantidAPI/BoxController.h"


constexpr size_t DIM_DISTRIBUTED_TEST = 3;


namespace Mantid {
namespace MDAlgorithms {

using BoxStructure = Mantid::DataObjects::MDBoxBase<Mantid::DataObjects::MDLeanEvent<DIM_DISTRIBUTED_TEST>, DIM_DISTRIBUTED_TEST>;


struct FrameInformation {
  std::string dimensionNames[3];
  Mantid::Kernel::SpecialCoordinateSystem specialCoordinateSystem;
  Mantid::Geometry::MDFrame_uptr frame;
};


struct BoxStructureInformation {
  std::unique_ptr<BoxStructure> boxStructure;
  Mantid::API::BoxController_sptr boxController;
};


class MANTID_MDALGORITHMS_DLL ConvertToDistributedMD : public API::ParallelAlgorithm {
public:
  using MDEventList = std::vector<Mantid::DataObjects::MDLeanEvent<DIM_DISTRIBUTED_TEST>>;

  const std::string name() const override;
  int version() const override;
  const std::string category() const override;
  const std::string summary() const override;

  static const std::string boxSplittingGroupName;

private:
  void init() override;
  void exec() override;

  BoxStructure* getPreliminaryBoxStructure(MDEventList& mdEvents) const;

  /// Send MDEvents to master
  void sendMDEventsToMaster(const Mantid::Parallel::Communicator& communicator,
                            const MDEventList& mdEvents) const;

  /// Receive MDEvents on master
  MDEventList receiveMDEventsOnMaster(const Mantid::Parallel::Communicator& communicator,
                               const MDEventList& mdEvents) const;

  BoxStructureInformation generatePreliminaryBoxStructure(const MDEventList& allEvents) const;

  void initBoxControllerProps(const std::string &SplitInto = "5",
                              int SplitThreshold = 1000,
                              int MaxRecursionDepth = 5);

  MDEventList getNPercentEvents(const Mantid::DataObjects::EventWorkspace& workspace) const;

  std::vector<coord_t> getWorkspaceExtents() const;

  BoxStructureInformation extractBoxStructure(Mantid::DataObjects::MDEventWorkspace3Lean& workspace) const;

  void sendRankResponsibility(const Mantid::Parallel::Communicator& communicator, const std::vector<std::pair<size_t, size_t>>& responsibility) const;

  std::vector<std::pair<size_t, size_t>> receiveRankResponsibility(const Mantid::Parallel::Communicator& communicator) const;

  // -------------------------------------------------------------------------------------------------------------------
  // Methods to build up a temporary MDEventWorkspace
  // -------------------------------------------------------------------------------------------------------------------
  FrameInformation createFrame() const;

  Mantid::DataObjects::MDEventWorkspace3Lean::sptr createTemporaryWorkspace() const;

  void addMDEventsToMDEventWorkspace(Mantid::DataObjects::MDEventWorkspace3Lean& workspace,
                                     const MDEventList & allEvents) const;


  void setBoxController(Mantid::API::BoxController_sptr bc) const;
  };

} // namespace MDAlgorithms
} // namespace Mantid

#endif /* MANTID_MDALGORITHMS_CONVERTTODISTRIBUTEDMD_H_ */
