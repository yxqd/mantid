#include "MantidAPI/BoxControllerAccess.h"
#include "MantidAPI/BoxController.h"

namespace Mantid {
namespace API {
void BoxControllerAccess::setNumMDBoxes(BoxController& boxController, const std::vector<size_t>& numMDBoxes) const {
  boxController.m_numMDBoxes = numMDBoxes;
}


void BoxControllerAccess::setNumMDGridBoxes(BoxController& boxController, const std::vector<size_t>& numMDGridBoxes) const {
  boxController.m_numMDGridBoxes = numMDGridBoxes;
}
} // namespace API
} // namespace Mantid
