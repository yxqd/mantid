#ifndef MANTID_MUON_MUONPAIRINGASYMMETRY_H_
#define MANTID_MUON_MUONPAIRINGASYMMETRY_H_

#include "MantidAPI/Algorithm.h"

using namespace Mantid::API;

namespace Mantid {
namespace Muon {

class DLLExport MuonPairingAsymmetry : public API::Algorithm {
public:
  MuonPairingAsymmetry() : API::Algorithm() {}
  ~MuonPairingAsymmetry() {}

  const std::string name() const override { return "MuonPairingAsymmetry"; }
  int version() const override { return (1); }
  const std::string category() const override { return "Muon\\DataHandling"; }
  const std::string summary() const override { return "."; }
  const std::vector<std::string> seeAlso() const override {
	  return{ "MuonProcess" };
  }

private:
  /// Initialisation code
  void init() override;
  /// Execution code
  void exec() override;
};

} // namespace Muon
} // namespace Mantid

#endif /* MANTID_MUON_MUONPAIRINGASYMMETRY_H_ */
