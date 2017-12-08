#ifndef MANTID_API_BOXCONTROLLERACCESS_H_
#define MANTID_API_BOXCONTROLLERACCESS_H_

#include "MantidAPI/DllConfig.h"
#include <vector>

namespace Mantid {
namespace API {

class BoxController;

class MANTID_API_DLL BoxControllerAccess {
public:
  void setNumMDBoxes(BoxController& boxController, const std::vector<size_t>& mdBoxes) const;
  void setNumMDGridBoxes(BoxController& boxController, const std::vector<size_t>& mdGridBoxes) const;
};

} // namespace API
} // namespace Mantid

#endif /* MANTID_API_BOXCONTROLLERACCESS_H_ */