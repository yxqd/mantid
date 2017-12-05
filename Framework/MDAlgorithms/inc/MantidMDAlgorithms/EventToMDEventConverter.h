#ifndef MANTID_MDALGORITHMS_EVENTTOMDEVENTCONVERTER_H_
#define MANTID_MDALGORITHMS_EVENTTOMDEVENTCONVERTER_H_

#include "MantidMDAlgorithms/DllConfig.h"
#include "MantidTypes/Core/DateAndTime.h"
#include "MantidTypes/Event/TofEvent.h"
#include "MantidDataObjects/MDEvent.h"
#include "MantidDataObjects/EventList.h"
#include "MantidAPI/SpectrumInfo.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidMDAlgorithms/MDWSDescription.h"

namespace Mantid {
namespace MDAlgorithms {

  constexpr size_t DIM_DISTRIBUTED_TEST = 3;


enum class QFrame { QSample, HKL, QLab };

struct EventConversionInfo {
  double l1 = 0.;
  double qSign = -1.;
  bool lorentzCorrection = false;
  Mantid::Kernel::V3D beamDirection;
  Mantid::Kernel::Matrix<double> ubMatrix;
  Mantid::Types::Core::DateAndTime cutOffTime;
  std::vector<coord_t> minExtents;
  std::vector<coord_t> maxExtents;
};

class MANTID_MDALGORITHMS_DLL EventToMDEventConverter {
public:
  std::vector<Mantid::DataObjects::MDLeanEvent<DIM_DISTRIBUTED_TEST>>
  getEvents(const Mantid::DataObjects::EventWorkspace &inputWorkspace,
            const std::vector<coord_t>& extents,
            double fraction = 1.,
            QFrame qFrame = QFrame::QLab,
            bool lorentzCorrection = false);

private:
  Mantid::Kernel::Matrix<double>
  getUbMatrix(const Mantid::DataObjects::EventWorkspace &inputWorkspace,
              QFrame qFrame) const;

  std::vector<Mantid::DataObjects::MDLeanEvent<DIM_DISTRIBUTED_TEST>>
  convertEvents(const size_t workspaceIndex,
                const Mantid::DataObjects::EventWorkspace &workspace,
                const EventConversionInfo &beamInformation,
                const API::SpectrumInfo &spectrumInfo);

  double getQSign() const;

  EventConversionInfo
  getEventConversionInfo(const Mantid::DataObjects::EventWorkspace &workspace,
                         QFrame qFrame, double fraction, bool lorentzCorrection,
                         const std::vector<coord_t>& extents);

  Mantid::Types::Core::DateAndTime
  getCutOffTime(const Mantid::DataObjects::EventWorkspace &workspace,
                double fraction) const;

  std::vector<Mantid::Types::Event::TofEvent>::const_iterator
  getEndIterator(const std::vector<Mantid::Types::Event::TofEvent> &events,
                 const EventConversionInfo &eventConversionInfo) const;

};

} // namespace MDAlgorithms
} // namespace Mantid

#endif /* MANTID_MDALGORITHMS_EVENTTOMDEVENTCONVERTER_H_ */
