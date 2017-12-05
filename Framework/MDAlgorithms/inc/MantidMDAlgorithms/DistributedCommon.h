#ifndef MANTID_MDALGORITHMS_DISTRIBUTEDCOMMON_H_
#define MANTID_MDALGORITHMS_DISTRIBUTEDCOMMON_H_

#include "MantidDataObjects/MDBoxBase.h"

namespace Mantid {
namespace MDAlgorithms {
namespace DistributedCommon {

constexpr size_t DIM_DISTRIBUTED_TEST = 3;
using BoxStructure = Mantid::DataObjects::MDBoxBase<
    Mantid::DataObjects::MDLeanEvent<DIM_DISTRIBUTED_TEST>,
    DIM_DISTRIBUTED_TEST>;
using MDEventList =
    std::vector<Mantid::DataObjects::MDLeanEvent<DIM_DISTRIBUTED_TEST>>;


struct BoxStructureInformation {
  std::unique_ptr<BoxStructure> boxStructure;
  Mantid::API::BoxController_sptr boxController;
};

}
}
}

#endif /* MANTID_MDALGORITHMS_BOXSTRUCTURESERIALIZER_H_ */