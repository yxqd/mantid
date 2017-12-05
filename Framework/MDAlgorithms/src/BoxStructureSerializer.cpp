#include "MantidMDAlgorithms/BoxStructureSerializer.h"
#include "MantidDataObjects/MDEventFactory.h"
#include <boost/make_shared.hpp>

namespace Mantid {
namespace MDAlgorithms {

using namespace DistributedCommon;

SerialBoxStructureInformation BoxStructureSerializer::serializeBoxStructure(
    const BoxStructureInformation &boxStructureInfo) const {

  SerialBoxStructureInformation serialBoxStructureInformation;
  initialize(serialBoxStructureInformation, boxStructureInfo.boxStructure.get(),
             *boxStructureInfo.boxController);
  return serialBoxStructureInformation;
}

BoxStructureInformation BoxStructureSerializer::deserializeBoxStructure(
    const SerialBoxStructureInformation &serialBoxStructureInformation) const {
  BoxStructureInformation boxStructureInformation;

  // Generate the new box controller
  const auto numberOfDimensions = serialBoxStructureInformation.numberOfDimensions;
  boxStructureInformation.boxController = boost::make_shared<Mantid::API::BoxController>(numberOfDimensions);
  boxStructureInformation.boxController->fromXMLString(serialBoxStructureInformation.boxController);
  auto boxController = boxStructureInformation.boxController;

  std::vector<API::IMDNode *> boxes;
  auto numberOfBoxes = serialBoxStructureInformation.boxType.size();
  boxes.assign(numberOfBoxes, nullptr);

  // We remove the check for the maximal number of dimensions and the event type.
  constexpr int eventType = 0;

  const auto& boxTypes = serialBoxStructureInformation.boxType;
  const auto& extents = serialBoxStructureInformation.extents;
  const auto& depth = serialBoxStructureInformation.depth;
  const auto& inverseVolume = serialBoxStructureInformation.inverseVolume;
  const auto& boxSignalErrorsquared = serialBoxStructureInformation.boxSignalErrorsquared;
  const auto& boxChildren = serialBoxStructureInformation.boxChildren;


  for (size_t i = 0; i < numberOfBoxes; i++) {
    API::IMDNode *box = nullptr;

    const auto boxType = boxTypes[i];
    if (boxType == 0)
      continue;

    // Extents
    std::vector<Mantid::Geometry::MDDimensionExtents<coord_t>> extentsVector(numberOfDimensions);
    for (size_t d = 0; d < numberOfDimensions; ++d) {
      extentsVector[d].setExtents(extents[i * numberOfDimensions * 2 + d * 2],
                                  extents[i * numberOfDimensions * 2 + d * 2 + 1]);
    }

    // Create the boxes
    if (boxType == 1) {
      box = DataObjects::MDEventFactory::createBox(numberOfDimensions, DataObjects::MDEventFactory::BoxType(eventType),
                                                   boxController, extentsVector, depth[i]);
    } else if (boxType == 2) {
      box = DataObjects::MDEventFactory::createBox(
        numberOfDimensions, DataObjects::MDEventFactory::BoxType(eventType+1), boxController,
        extentsVector, depth[i]);
    } else {
      continue;
    }

    // Set the id
    box->setID(i);
    box->calcVolume();

    auto volume = inverseVolume[i];
    if (volume <= FLT_EPSILON) {
      volume = 1;
    }

    if (std::fabs((box->getInverseVolume() - volume) / volume) > 1.e-5) {
      box->setInverseVolume(static_cast<coord_t>(inverseVolume[i]));
    }

    // Set the cached values
    box->setSignal(boxSignalErrorsquared[i * 2]);
    box->setErrorSquared(boxSignalErrorsquared[i * 2 + 1]);

    // Save the box at its index in the vector.
    boxes[i] = box;
  }

  // Set the parents fo the children
  for (size_t i = 0; i < numberOfBoxes; i++) {
    if (boxTypes[i] == 2) {
      size_t indexStart = static_cast<size_t>(boxChildren[i * 2]);
      size_t indexEnd = static_cast<size_t>(boxChildren[i * 2 + 1] + 1);
      boxes[i]->setChildren(boxes, indexStart, indexEnd);
    }
  }
  boxController->setMaxId(numberOfBoxes);

  // Set everything on the
  auto rootNode = dynamic_cast<BoxStructure*>(boxes[0]);
  if (!rootNode) {
    throw std::runtime_error("The conversion from INode to MDBoxBase was not successful");
  }

  std::unique_ptr<BoxStructure> uniqueBoxStructure(rootNode);
  boxStructureInformation.boxStructure = std::move(uniqueBoxStructure);
  return boxStructureInformation;
}

void BoxStructureSerializer::initialize(
    SerialBoxStructureInformation &serialBoxStructureInformation,
    BoxStructure *boxStructure,
    const Mantid::API::BoxController &boxController) const {
  // Set the box controller
  serialBoxStructureInformation.boxController = boxController.toXMLString();
  const auto numberOfDimensions = boxController.getNDims();
  serialBoxStructureInformation.numberOfDimensions = numberOfDimensions;
  // Initialize the vectors
  std::vector<Mantid::API::IMDNode *> boxes;
  boxStructure->getBoxes(boxes, 1000, false);
  API::IMDNode::sortObjByID(boxes);

  // Reserve size
  auto maxBoxes = boxes.size();
  serialBoxStructureInformation.boxType.assign(maxBoxes, 0);
  serialBoxStructureInformation.depth.assign(maxBoxes, -1);
  serialBoxStructureInformation.boxEventIndex.assign(maxBoxes * 2, 0);
  serialBoxStructureInformation.extents.assign(
      maxBoxes * numberOfDimensions * 2, 0);
  serialBoxStructureInformation.inverseVolume.assign(maxBoxes, 0);
  serialBoxStructureInformation.boxSignalErrorsquared.assign(maxBoxes * 2, 0);
  serialBoxStructureInformation.boxChildren.assign(maxBoxes * 2, 0);

  // Populate the serialization vector
  API::IMDNode *box;
  size_t ic(0);
  bool filePositionDefined(true);
  for (size_t i = 0; i < maxBoxes; i++) {
    box = boxes[i];
    auto id = box->getID();
    auto numChildren = box->getNumChildren();
    if (numChildren > 0) {
      serialBoxStructureInformation.boxType[ic] = 2;
      serialBoxStructureInformation.boxChildren[ic * 2] =
          int(box->getChild(0)->getID());
      serialBoxStructureInformation.boxChildren[ic * 2 + 1] =
          int(box->getChild(numChildren - 1)->getID());
      serialBoxStructureInformation.boxEventIndex[ic * 2] = 0;
      serialBoxStructureInformation.boxEventIndex[ic * 2 + 1] = 0;
    } else {
      serialBoxStructureInformation.boxType[ic] = 1;
      serialBoxStructureInformation.boxChildren[ic * 2] = 0;
      serialBoxStructureInformation.boxChildren[ic * 2 + 1] = 0;
      auto nPoints = box->getNPoints();
      Kernel::ISaveable *pSaver = box->getISaveable();
      if (pSaver)
        serialBoxStructureInformation.boxEventIndex[ic * 2] =
            pSaver->getFilePosition();
      else
        filePositionDefined = false;

      serialBoxStructureInformation.boxEventIndex[ic * 2 + 1] = nPoints;
    }

    serialBoxStructureInformation.depth[ic] = int(box->getDepth());
    serialBoxStructureInformation.boxSignalErrorsquared[ic * 2] =
        double(box->getSignal());
    serialBoxStructureInformation.boxSignalErrorsquared[ic * 2 + 1] =
        double(box->getErrorSquared());
    serialBoxStructureInformation.inverseVolume[ic] = box->getInverseVolume();
    for (size_t d = 0; d < numberOfDimensions; d++) {
      size_t newIndex = id * size_t(numberOfDimensions * 2) + d * 2;
      serialBoxStructureInformation.extents[newIndex] =
          box->getExtents(d).getMin();
      serialBoxStructureInformation.extents[newIndex + 1] =
          box->getExtents(d).getMax();
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
