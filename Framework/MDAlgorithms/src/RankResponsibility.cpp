#include <MantidAPI/IMDNode.h>
#include "MantidMDAlgorithms/RankResponsibility.h"


namespace Mantid {
namespace MDAlgorithms {

std::vector<std::pair<size_t, size_t>> RankResponsibility::getResponsibilites(int numberOfRanks, BoxStructure* boxStructure) const {
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
  for (size_t index = 0; index < boxes.size(); ++index) {
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

    // Count the data
    signalOnCurrentRank += boxes[index]->getSignal();
  }

  // The last rank has not been written out
  boxResponsibilityRangeOnRank.emplace_back(startIndex, boxes.size()-1);

  return boxResponsibilityRangeOnRank;
}


} // namespace MDAlgorithms
} // namespace Mantid
