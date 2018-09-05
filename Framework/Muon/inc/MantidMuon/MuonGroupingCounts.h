#ifndef MANTID_MUON_MUONGROUPINGCOUNTS_H_
#define MANTID_MUON_MUONGROUPINGCOUNTS_H_

#include "MantidAPI/Algorithm.h"

using namespace Mantid::API;

namespace Mantid {
namespace Muon {

class DLLExport MuonGroupingCounts : public API::Algorithm {
public:
  MuonGroupingCounts() : API::Algorithm() {}
  ~MuonGroupingCounts() {}

  const std::string name() const override { return "MuonGroupingCounts"; }
  int version() const override { return (1); }
  const std::string category() const override { return "Muon\\DataHandling"; }
  const std::string summary() const override { return "."; }
  const std::vector<std::string> seeAlso() const override {
    return {"MuonProcess"};
  }

private:
  /// Initialisation code
  void init() override;
  /// Execution code
  void exec() override;
};

} // namespace Muon
} // namespace Mantid

#endif /* MANTID_MUON_MUONGROUPINGCOUNTS_H_ */
