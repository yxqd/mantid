#include <MantidAPI/IMDNode.h>
#include "MantidMDAlgorithms/RankResponsibility.h"


namespace Mantid {
namespace MDAlgorithms {

std::vector<std::pair<size_t, size_t>> RankResponsibility::getResponsibilites(int numberOfRanks, BoxStructure* boxStructure) const {
  // The ideal partition strategy is able to split the number of events, ie the signal, evenly between the ranks.
  // However the signal is contain in boxes which we cannot partition further, this means that the minimal unit
  // of events is one box.

  // TODO: After a first pass of the partitioning we should check if it is an even an fair split. If this is not the
  //       case then we should modify the split parameters, e.g. lower the split threshold and so on. This is
  //       2nd order work though.

  // We determine how many events should be placed into each rank
  const auto totalSignal = boxStructure->getSignal();
  const auto signalPerRank = totalSignal/numberOfRanks;


  // We get only the leaf nodes
  std::vector<Mantid::API::IMDNode*> boxes;
  boxStructure->getBoxes(boxes, 1000, true /* only leaf nodes*/);

  size_t startIndex = 0;
  signal_t signalOnCurrentRank = 0.;
  std::vector<std::pair<size_t, size_t>> boxResponsibilityRangeOnRank;
  int currentRank = 0;

  // The approach is to iterate through the boxes and sum up the signal until we exceed our threshold (signalPerRank).
  // We then reset our sum counter and associate the next range with the next rank and so on.
  // We need to ensure though that we do not end up with ranks without any boxes. This can happen:
  // i. Because the load is highly concentrated in very few boxes and there are many ranks.
  // ii. There are more ranks than boxes (in this case we throw an exception for now, in the future we should resplit)

  if (numberOfRanks > boxes.size()) {
    throw std::runtime_error("There are more ranks and boxes. This cannot be handled currently");
  }

  auto numberOfBoxes = boxes.size();
  for (size_t index = 0; index < numberOfBoxes; ++index) {
    auto numberUnassignedBoxes = numberOfBoxes - index;
    auto numberUnassignedRanks = numberOfRanks - currentRank;

    // If we have an equal number of unassigned ranks and boxes, then assign one box per rank
    if (numberUnassignedBoxes == static_cast<size_t>(numberUnassignedRanks)) {
      if (currentRank != numberOfRanks - 1) {
        boxResponsibilityRangeOnRank.emplace_back(startIndex, index-1);
        startIndex = index;
        signalOnCurrentRank = 0;
        ++currentRank;
      }
    } else {
      // Perform standard assignment
      // If we have reached the max amount of signal per rank,
      // then distribute to the next rank
      if (signalOnCurrentRank >= signalPerRank) {
        // If we are already at the last rank then we leave it there for now
        if (currentRank != numberOfRanks - 1) {
          boxResponsibilityRangeOnRank.emplace_back(startIndex, index-1);
          startIndex = index;
          signalOnCurrentRank = 0;
          ++currentRank;
        }
      }
    }
    // Count the data
    signalOnCurrentRank += boxes[index]->getSignal();
  }

  // The last rank has not been written out
  boxResponsibilityRangeOnRank.emplace_back(startIndex, boxes.size()-1);

  return boxResponsibilityRangeOnRank;
}


} // namespace MDAlgorithms
} // namespace Mantid
