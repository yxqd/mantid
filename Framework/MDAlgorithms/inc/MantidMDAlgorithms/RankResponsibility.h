#ifndef MANTID_MDALGORITHMS_RANKRESPONSIBILITY_H_
#define MANTID_MDALGORITHMS_RANKRESPONSIBILITY_H_

#include "MantidMDAlgorithms/DllConfig.h"
#include "MantidDataObjects/MDLeanEvent.h"
#include "MantidDataObjects/MDBoxBase.h"

#include <vector>

namespace Mantid {
namespace MDAlgorithms {


class MANTID_MDALGORITHMS_DLL RankResponsibility {
public:
  using BoxStructure = Mantid::DataObjects::MDBoxBase<Mantid::DataObjects::MDLeanEvent<3>, 3>;

  std::vector<std::pair<size_t, size_t>> getResponsibilites(int numberOfRanks, BoxStructure* boxStructure) const;

};

} // namespace MDAlgorithms
} // namespace Mantid

#endif /* MANTID_MDALGORITHMS_RANKRESPONSIBILITY_H_ */