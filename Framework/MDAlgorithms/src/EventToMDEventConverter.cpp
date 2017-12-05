#include "MantidMDAlgorithms/EventToMDEventConverter.h"
#include "MantidAPI/AlgorithmManager.h"
#include "MantidAPI/Run.h"
#include "MantidAPI/Sample.h"
#include "MantidKernel/ConfigService.h"

using namespace Mantid::Kernel;
using namespace Mantid::DataObjects;
using namespace Mantid::Types::Event;
using namespace Mantid::Types::Core;

namespace Mantid {
namespace MDAlgorithms {

using Lev3D = Mantid::DataObjects::MDLeanEvent<DIM_DISTRIBUTED_TEST>;

std::vector<Lev3D>
EventToMDEventConverter::getEvents(const EventWorkspace &workspace,
                                   const std::vector<coord_t>& extents,
                                   double fraction,
                                   QFrame qFrame,
                                   bool lorentzCorrection) {
  // Get beamline inforamtion
  auto eventConversionInfo =
      getEventConversionInfo(workspace, qFrame, fraction, lorentzCorrection, extents);

  // Convert the number of events
  const auto &spectrumInfo = workspace.spectrumInfo();
  const auto totalNumberOfSpectra = workspace.getNumberHistograms();
  const auto totalNumberOfEvents = workspace.getNumberEvents();
  std::vector<Lev3D> output;
  output.reserve(totalNumberOfEvents);
  for (auto workspaceIndex = 0ul; workspaceIndex < totalNumberOfSpectra;
       ++workspaceIndex) {
    // Get Event list
    auto events = convertEvents(workspaceIndex, workspace, eventConversionInfo,
                                spectrumInfo);
    std::move(events.begin(), events.end(), std::back_inserter(output));
  }

  // Maybe we need to shrink (if we filter out events for example)
  return output;
}

Mantid::Kernel::Matrix<double>
EventToMDEventConverter::getUbMatrix(const EventWorkspace &workspace,
                                     QFrame qFrame) const {
  auto ubMatrix = Matrix<double>(3, 3, true);

  switch (qFrame) {
  case QFrame::QSample: {
    ubMatrix = workspace.run().getGoniometerMatrix();
    ubMatrix.Invert();
    break;
  }
  case QFrame::HKL: {
    const auto &ub = workspace.sample().getOrientedLattice().getUB();
    const auto &goniometer = workspace.run().getGoniometerMatrix();
    ubMatrix = goniometer * ub;
    ubMatrix.Invert();
    ubMatrix /= (2 * M_PI);
    break;
  }
  case QFrame::QLab: {
    break;
  }
  }
  return ubMatrix;
}

std::vector<Lev3D> EventToMDEventConverter::convertEvents(
    const size_t workspaceIndex, const EventWorkspace &workspace,
    const EventConversionInfo &eventConversionInfo,
    const API::SpectrumInfo &spectrumInfo) {
  const EventList &eventList = workspace.getSpectrum(workspaceIndex);
  std::vector<Lev3D> mdEvents;
  // Get the position of the detector there.
  const auto &detectors = eventList.getDetectorIDs();
  if (!detectors.empty()) {
    // Check if a detector is located at this workspace index, returns
    // immediately if one is not found.
    if (!spectrumInfo.hasDetectors(workspaceIndex)) {
      return std::vector<Lev3D>();
    }

    // Neutron's total travelled distance
    const auto distance =
        eventConversionInfo.l1 + spectrumInfo.l2(workspaceIndex);

    // Vector between the sample and the detector
    const auto detPos = spectrumInfo.position(workspaceIndex);

    // Detector direction normalized to 1
    const auto detDir = detPos / detPos.norm();

    // The direction of momentum transfer in the inelastic convention ki-kf
    //  = input beam direction (normalized to 1) - output beam direction
    //  (normalized to 1)
    auto qDirLabFrame = eventConversionInfo.beamDirection - detDir;

    // Take into account the crystallographic convention
    qDirLabFrame *= eventConversionInfo.qSign;

    // Convert into the right frame
    auto qDir = eventConversionInfo.ubMatrix * qDirLabFrame;

    // For speed we extract the components.
    coord_t Q_dir_x = coord_t(qDir.X());
    coord_t Q_dir_y = coord_t(qDir.Y());
    coord_t Q_dir_z = coord_t(qDir.Z());

    // For lorentz correction, calculate  sin(theta))^2
    auto sinThetaSquared = 0.;
    if (eventConversionInfo.lorentzCorrection) {
      // Scattering angle = 2 theta = angle between neutron beam direction and
      // the detector (scattering) direction
      // The formula for Lorentz Correction is sin(theta), i.e. sin(half the
      // scattering angle)
      const auto theta = spectrumInfo.twoTheta(workspaceIndex) / 2.0;
      sinThetaSquared = sin(theta);
      sinThetaSquared = sinThetaSquared * sinThetaSquared;
    }

    /** Constant that you divide by tof (in usec) to get wavenumber in ang^-1 :
     * Wavenumber (in ang^-1) =  (PhysicalConstants::NeutronMass * distance) /
     * ((tof (in usec) * 1e-6) * PhysicalConstants::h_bar) * 1e-10; */
    const auto conversionFactor =
        (PhysicalConstants::NeutronMass * distance * 1e-10) /
        (1e-6 * PhysicalConstants::h_bar);

    // This little dance makes the getting vector of events more general (since
    // you can't overload by return type).
    typename std::vector<Mantid::Types::Event::TofEvent> const *events_ptr;
    getEventsFrom(eventList, events_ptr);
    typename std::vector<Mantid::Types::Event::TofEvent> const &events =
        *events_ptr;

    // Iterators to start/end
    auto it = events.begin();
    auto it_end = getEndIterator(events, eventConversionInfo);

    const auto numberOfEvents = eventList.getNumberEvents();
    mdEvents.reserve(numberOfEvents);

    for (; it != it_end; it++) {
      // Get the wavenumber in ang^-1 using the previously calculated constant.
      coord_t wavenumber = coord_t(conversionFactor / it->tof());

      // Q vector = K_final - K_initial = wavenumber * (output_direction -
      // input_direction)
      coord_t center[3] = {Q_dir_x * wavenumber, Q_dir_y * wavenumber,
                           Q_dir_z * wavenumber};

//      // Check that the event is within bounds
//      if (center[0] < eventConversionInfo.minExtents[0] || center[0] >= eventConversionInfo.maxExtents[0])
//        continue;
//      if (center[1] < eventConversionInfo.minExtents[1] || center[1] >=  eventConversionInfo.maxExtents[1])
//        continue;
//      if (center[2] < eventConversionInfo.minExtents[2] || center[2] >= eventConversionInfo.maxExtents[2])
//        continue;

      if (eventConversionInfo.lorentzCorrection) {
        // double lambda = 1.0/wavenumber;
        // (sin(theta))^2 / wavelength^4
        auto correct = float(sinThetaSquared * wavenumber * wavenumber *
                              wavenumber * wavenumber);
        // Push the MDLeanEvent but correct the weight.
        mdEvents.emplace_back(float(it->weight() * correct),
                              float(it->errorSquared() * correct * correct),
                              center);
      } else {
        // Push the MDLeanEvent with the same weight
        mdEvents.emplace_back(float(it->weight()), float(it->errorSquared()),
                              center);
      }
    }
  }
  return mdEvents;
}

std::vector<TofEvent>::const_iterator EventToMDEventConverter::getEndIterator(
    const std::vector<TofEvent> &events,
    const EventConversionInfo &eventConversionInfo) const {
  if (events.empty()) {
    return events.cend();
  }

  if (events.back().pulseTime() <= eventConversionInfo.cutOffTime) {
    return events.cend();
  } else {
    // We need to find the first element in the events list which has a pulse
    // time which is larger than the cutoff time
    const auto &cutOffTime = eventConversionInfo.cutOffTime;
    auto comparePulseTime =
        [](const TofEvent &event, const DateAndTime &cutOff)
            -> bool { return event.pulseTime() < cutOff; };

    const auto largerOrEqual = std::lower_bound(events.begin(), events.end(),
                                          cutOffTime, comparePulseTime);
    return largerOrEqual != events.cend() &&
                   !comparePulseTime(*largerOrEqual, cutOffTime)
               ? largerOrEqual
               : events.cend();
  }
}

double EventToMDEventConverter::getQSign() const {
  auto qSign = -1.0;
  const auto convention =
      Mantid::Kernel::ConfigService::Instance().getString("Q.convention");
  if (convention == "Crystallography")
    qSign = 1.0;
  return qSign;
}

EventConversionInfo EventToMDEventConverter::getEventConversionInfo(
    const Mantid::DataObjects::EventWorkspace &workspace, QFrame qFrame,
    double fraction, bool lorentzCorrection, const std::vector<coord_t>& extents) {
  const auto &spectrumInfo = workspace.spectrumInfo();
  EventConversionInfo eventConversionInfo;
  eventConversionInfo.lorentzCorrection = lorentzCorrection;

  // Get general beam information
  eventConversionInfo.l1 = spectrumInfo.l1();
  const auto samplePosition = spectrumInfo.samplePosition();
  const auto sourcePosition = spectrumInfo.sourcePosition();
  const auto beamline = samplePosition - sourcePosition;
  eventConversionInfo.beamDirection = beamline / beamline.norm();

  // Get the UB matrix
  eventConversionInfo.ubMatrix = getUbMatrix(workspace, qFrame);

  // Set the QSign
  eventConversionInfo.qSign = getQSign();

  // Get the start and stop time
  auto cutOffTime = getCutOffTime(workspace, fraction);
  eventConversionInfo.cutOffTime = cutOffTime;

  // Populate the extents
  if (extents.size() != 2*DIM_DISTRIBUTED_TEST) {
    throw std::runtime_error("The number extents has to be twice the number of dimensions");
  }

  for(auto dimension=0ul; dimension < DIM_DISTRIBUTED_TEST; ++dimension) {
    eventConversionInfo.minExtents.emplace_back(extents[dimension]);
    eventConversionInfo.maxExtents.emplace_back(extents[dimension+1]);
  }

  return eventConversionInfo;
}

Mantid::Types::Core::DateAndTime EventToMDEventConverter::getCutOffTime(
    const Mantid::DataObjects::EventWorkspace &workspace,
    double fraction) const {
  auto maxTime = workspace.getPulseTimeMax();

  if (fraction == 1.f) {
    return maxTime;
  }

  auto minTime = workspace.getPulseTimeMin();
  auto duration = maxTime - minTime;
  auto durationNano = DateAndTime::nanosecondsFromDuration(duration);
  auto reducedDurationNanoSeconds =
      static_cast<int64_t>(static_cast<float>(durationNano) * fraction);
  return minTime + reducedDurationNanoSeconds;
}

} // namespace MDAlgorithms
} // namespace Mantid
