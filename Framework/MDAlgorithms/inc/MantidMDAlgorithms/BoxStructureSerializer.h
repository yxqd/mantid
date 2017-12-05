#ifndef MANTID_MDALGORITHMS_BOXSTRUCTURESERIALIZER_H_
#define MANTID_MDALGORITHMS_BOXSTRUCTURESERIALIZER_H_

#include "MantidMDAlgorithms/DllConfig.h"
#include "MantidDataObjects/MDLeanEvent.h"
#include "MantidAPI/BoxController.h"
#include "MantidDataObjects/MDBoxBase.h"

#include <boost/serialization/access.hpp>

#include <vector>
#include <string>

namespace Mantid {
namespace MDAlgorithms {


struct SerialBoxStructureInformation {
private:
  friend class boost::serialization::access;
public:
  template<class Archive>
  void serialize(Archive &ar, const unsigned int /*version*/)
  {
    ar & boxType;
    ar & depth;
    ar & boxEventIndex;
    ar & extents;
    ar & inverseVolume;
    ar & boxSignalErrorsquared;
    ar & boxChildren;
    ar & numberOfDimensions;
    ar & boxController;
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
  using BoxStructure = Mantid::DataObjects::MDBoxBase<Mantid::DataObjects::MDLeanEvent<3>, 3>;

  SerialBoxStructureInformation createSerializedBoxStructure(BoxStructure* boxStructure, const Mantid::API::BoxController& boxController
  ) const;

private:
  void initialize(SerialBoxStructureInformation& serialBoxStructureInformation, BoxStructure* boxStructure, const Mantid::API::BoxController& boxController) const;

};

} // namespace MDAlgorithms
} // namespace Mantid

#endif /* MANTID_MDALGORITHMS_BOXSTRUCTURESERIALIZER_H_ */