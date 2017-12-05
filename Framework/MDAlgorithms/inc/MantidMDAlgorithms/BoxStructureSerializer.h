#ifndef MANTID_MDALGORITHMS_BOXSTRUCTURESERIALIZER_H_
#define MANTID_MDALGORITHMS_BOXSTRUCTURESERIALIZER_H_

#include "MantidMDAlgorithms/DllConfig.h"
#include "MantidMDAlgorithms/DistributedCommon.h"
#include "MantidAPI/BoxController.h"
#include "MantidDataObjects/MDLeanEvent.h"
#include "MantidDataObjects/MDBoxBase.h"

#include <boost/serialization/access.hpp>

#include <string>


namespace Mantid {
namespace MDAlgorithms {

struct SerialBoxStructureInformation {
private:
  friend class boost::serialization::access;

public:
  template <class Archive>
  void serialize(Archive &ar, const unsigned int /*version*/) {
    ar &boxType;
    ar &depth;
    ar &boxEventIndex;
    ar &extents;
    ar &inverseVolume;
    ar &boxSignalErrorsquared;
    ar &boxChildren;
    ar &numberOfDimensions;
    ar &boxController;
  }

  std::vector<int> boxType;
  std::vector<int> depth;
  std::vector<uint64_t> boxEventIndex;
  std::vector<double> extents;
  std::vector<double> inverseVolume;
  std::vector<double> boxSignalErrorsquared;
  std::vector<int> boxChildren;
  size_t numberOfDimensions;
  std::string boxController;
};

class MANTID_MDALGORITHMS_DLL BoxStructureSerializer {
public:
  SerialBoxStructureInformation serializeBoxStructure(const DistributedCommon::BoxStructureInformation& boxStructureInfo) const;

  DistributedCommon::BoxStructureInformation deserializeBoxStructure(const SerialBoxStructureInformation& serialBoxStructureInformation) const;

private:
  void initialize(SerialBoxStructureInformation &serialBoxStructureInformation,
                  DistributedCommon::BoxStructure *boxStructure,
                  const Mantid::API::BoxController &boxController) const;
};

} // namespace MDAlgorithms
} // namespace Mantid

#endif /* MANTID_MDALGORITHMS_BOXSTRUCTURESERIALIZER_H_ */