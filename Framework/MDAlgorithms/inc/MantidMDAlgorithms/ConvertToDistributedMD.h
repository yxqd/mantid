#ifndef MANTID_MDALGORITHMS_CONVERTTODISTRIBUTEDMD_H_
#define MANTID_MDALGORITHMS_CONVERTTODISTRIBUTEDMD_H_

#include "MantidMDAlgorithms/DllConfig.h"
#include "MantidAPI/ParallelAlgorithm.h"
#include "MantidDataObjects/MDEvent.h"
#include "MantidDataObjects/MDBoxBase.h"
#include "MantidDataObjects/MDLeanEvent.h"
#include "MantidAPI/ParallelAlgorithm.h"


constexpr size_t DIM_DISTRIBUTED_TEST = 3;


namespace Mantid {
namespace MDAlgorithms {

class MANTID_MDALGORITHMS_DLL ConvertToDistributedMD : public API::ParallelAlgorithm {
public:
  using BoxStructure = Mantid::DataObjects::MDBoxBase<Mantid::DataObjects::MDLeanEvent<DIM_DISTRIBUTED_TEST>, DIM_DISTRIBUTED_TEST>;
  using MDEventList = std::vector<Mantid::DataObjects::MDLeanEvent<DIM_DISTRIBUTED_TEST>>;

  const std::string name() const override;
  int version() const override;
  const std::string category() const override;
  const std::string summary() const override;

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

//  BoxStructure* generatePreliminaryBoxStructure(const MDEventList& allEvents) const;
};

} // namespace MDAlgorithms
} // namespace Mantid

#endif /* MANTID_MDALGORITHMS_CONVERTTODISTRIBUTEDMD_H_ */
