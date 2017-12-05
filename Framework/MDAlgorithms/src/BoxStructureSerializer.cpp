#include "MantidMDAlgorithms/BoxStructureSerializer.h"

namespace Mantid {
namespace MDAlgorithms {

SerialBoxStructureInformation BoxStructureSerializer::createSerializedBoxStructure(BoxStructure* boxStructure,
                                                                                   const Mantid::API::BoxController& boxController) const {

  SerialBoxStructureInformation serialBoxStructureInformation;
  initialize(serialBoxStructureInformation, boxStructure, boxController);
  return serialBoxStructureInformation;
}

void BoxStructureSerializer::initialize(SerialBoxStructureInformation& serialBoxStructureInformation,
                                        BoxStructure* boxStructure,
                                        const Mantid::API::BoxController& boxController) const {
  // Set the box controller
  serialBoxStructureInformation.boxController = boxController.toXMLString();
  const auto numberOfDimensions = boxController.getNDims();
  serialBoxStructureInformation.numberOfDimensions = numberOfDimensions;
  // Initialize the vectors
  std::vector<Mantid::API::IMDNode*> boxes;
  boxStructure->getBoxes(boxes, 1000, false);
  API::IMDNode::sortObjByID(boxes);

  // Reserve size
  auto maxBoxes = boxes.size();
  serialBoxStructureInformation.boxType.assign(maxBoxes, 0);
  serialBoxStructureInformation.depth.assign(maxBoxes, -1);
  serialBoxStructureInformation.boxEventIndex.assign(maxBoxes*2, 0);
  serialBoxStructureInformation.extents.assign(maxBoxes*numberOfDimensions*2, 0);
  serialBoxStructureInformation.inverseVolume.assign(maxBoxes, 0);
  serialBoxStructureInformation.boxSignalErrorsquared.assign(maxBoxes*2, 0);
  serialBoxStructureInformation.boxChildren.assign(maxBoxes*2, 0);

  // Populate the serialization vector
  API::IMDNode *box;
  size_t ic(0);
  bool filePositionDefined(true);
  for (size_t i = 0; i < maxBoxes; i++) {
    box = boxes[i];
    auto id = box->getID();
    auto numChildren = box->getNumChildren();
    if (numChildren > 0)
    {
      serialBoxStructureInformation.boxType[ic] = 2;
      serialBoxStructureInformation.boxChildren[ic * 2] = int(box->getChild(0)->getID());
      serialBoxStructureInformation.boxChildren[ic * 2 + 1] = int(box->getChild(numChildren - 1)->getID());
      serialBoxStructureInformation.boxEventIndex[ic * 2] = 0;
      serialBoxStructureInformation.boxEventIndex[ic * 2 + 1] = 0;
    } else {
      serialBoxStructureInformation.boxType[ic] = 1;
      serialBoxStructureInformation.boxChildren[ic * 2] = 0;
      serialBoxStructureInformation.boxChildren[ic * 2 + 1] = 0;
      auto nPoints = box->getNPoints();
      Kernel::ISaveable *pSaver = box->getISaveable();
      if (pSaver)
        serialBoxStructureInformation.boxEventIndex[ic * 2] = pSaver->getFilePosition();
      else
        filePositionDefined = false;

      serialBoxStructureInformation.boxEventIndex[ic * 2 + 1] = nPoints;
    }

    serialBoxStructureInformation.depth[ic] = int(box->getDepth());
    serialBoxStructureInformation.boxSignalErrorsquared[ic * 2] = double(box->getSignal());
    serialBoxStructureInformation.boxSignalErrorsquared[ic * 2 + 1] = double(box->getErrorSquared());
    serialBoxStructureInformation.inverseVolume[ic] = box->getInverseVolume();
    for (size_t d = 0; d < numberOfDimensions; d++) {
      size_t newIndex = id * size_t(numberOfDimensions * 2) + d * 2;
      serialBoxStructureInformation.extents[newIndex] = box->getExtents(d).getMin();
      serialBoxStructureInformation.extents[newIndex + 1] = box->getExtents(d).getMax();
    }
    ic++;
  }
  // file position have to be calculated afresh
  if (!filePositionDefined) {
    uint64_t boxPosition(0);
    for (size_t i = 0; i < maxBoxes; i++) {
      if (serialBoxStructureInformation.boxType[i] == 1) {
        serialBoxStructureInformation.boxEventIndex[2 * i] = boxPosition;
        boxPosition += serialBoxStructureInformation.boxEventIndex[2 * i + 1];
      }
    }
  }
}


} // namespace MDAlgorithms
} // namespace Mantid
