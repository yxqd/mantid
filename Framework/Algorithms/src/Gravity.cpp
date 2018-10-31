#include "MantidAlgorithms/Gravity.h"
#include "MantidAPI/FrameworkManager.h"
#include "MantidAPI/AnalysisDataService.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/TextAxis.h"
#include "MantidAPI/NumericAxis.h"
#include "MantidKernel/MersenneTwister.h"
#include "MantidKernel/PseudoRandomNumberGenerator.h"
#include "MantidKernel/normal_distribution.h"
#include "MantidGeometry/Instrument.h"
#include "MantidGeometry/Instrument/DetectorInfo.h"
#include "MantidGeometry/Instrument/CompAssembly.h"
#include <Eigen/LU>

#include <iostream>

using namespace Mantid::API;
using namespace Mantid::Geometry;
using namespace Mantid::Kernel;

namespace Mantid {
namespace Algorithms {

DECLARE_ALGORITHM(Gravity)

const std::string Gravity::name() const { return "Gravity"; }

int Gravity::version() const { return 1; }

const std::string Gravity::category() const { return "Arithmetic"; }

const std::string Gravity::summary() const { return "Computes gravity."; }

//----------------------------------------------------------------------------------------------
namespace {

double const h_plank = 6.63e-34;
double const m_neutron = 1.67e-27;
double const lam_to_speed = h_plank / m_neutron * 1e10;
double const G = 9.80665;
// random number generator
std::mt19937 rng;

class TrajectoryError: public std::runtime_error {
public:
  TrajectoryError(std::string const&msg) : std::runtime_error(msg) {}
};

//----------------------------------------------------------------------------------------------
/// Fill in a spectrum of a MatrixWorkspace
template <typename Generator>
void fillSpectrum(MatrixWorkspace &ws, size_t index, double t0, double t1, Generator gen) {
    auto &x = ws.mutableX(index);
    auto &y = ws.mutableY(index);
    auto nBins = ws.blocksize();
    auto dt = (t1 - t0) / (nBins - 1);
    for(size_t i = 0; i < nBins; ++i) {
      auto t = t0 + double(i) * dt;
      x[i] = gen.x(t);
      y[i] = gen.y(t);
    }
}

struct Plane {
  double xm, ym, theta;
  double x(double t) const { return t * cos(theta) + xm; }
  double y(double t) const { return t * sin(theta) + ym; }
};

struct Vertical {
  double xm;
  double x(double t) const { return xm; }
  double y(double t) const { return t; }
};

//----------------------------------------------------------------------------------------------

/// Calculate time that a particle takes to hit an infinte line. The line passes through point (0, 0) and
///    inclined to the horizon at some angle.
///    Coordinates: x - horizontal, y - vertical, up.
/// @param vx :: x component of the initial velocity.
/// @param vy :: y component of the initial velocity.
/// @param g :: Gravitational constant.
/// @param c :: Cosine of the angle of inclination of the line.
/// @param s :: Sine of the angle of inclination of the line.
/// @param x0 :: Initial x coordinate of the particle.
/// @param y0 :: Initial y coordinate of the particle.
double calc_t1(double vx, double vy, double g, double c, double s, double x0, double y0, double xm , double ym) {
  double const d = 2*c*c*g*y0 - 2*c*c*g*ym + c*c*vy*vy - 2*c*g*s*x0 + 2*c*g*s*xm - 2*c*s*vx*vy + s*s*vx*vx;
  if (d < 0) {
    throw TrajectoryError("Particle cannot hit the plane.");
  }
  double const r = c*vy - s*vx;
  auto const t1 = (r - sqrt(d))/(c*g);
  auto const t2 = (r + sqrt(d))/(c*g);
  if (t2 > 0 && t1 > 0) {
    return t1 < t2 ? t1 : t2;
  } else if (t1 > 0) {
    return t1;
  } else if (t2 > 0) {
    return t2;
  }
  throw TrajectoryError("Particle cannot hit the plane.");
}

/// A path of a free falling particle
/// x0, y0 - initial coordinates
/// vx, vy - initial velocity
/// x1, y1 - final coordinates
/// v1x, v1y - final velocity
/// t - time it takes to get from the initial to the final point
/// Invariants:
///    path.x1 == path.x(path.t)
///    path.y1 == path.y(path.t)
///    path.v1x == path.vx
///    path.v1y == path.vy - G * path.t
struct Path {
  double t, x0, y0, vx, vy, x1, y1, v1x, v1y;
  double x(double tt) const { return x0 + vx * tt; }
  double y(double tt) const { return y0 + tt * (vy - G * tt / 2.0); }
  double displacementLength2() const {return (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0);}
};

class Trajectory: public std::vector<Path> {
public:
  Trajectory() = default;
  Trajectory(std::initializer_list<Path> l) : std::vector<Path>(l) {}
};

double calcAngleToPassBetween2Points(double x0, double y0, double x1, double y1, double v) {
  auto const dx = x1 - x0;
  auto const dy = y1 - y0;
  auto const v2 = v * v;
  auto const D = -G*G*dx*dx + v2*(-2*G*dy + v2);
  if (D < 0.0) {
    throw TrajectoryError("Particle cannot pass through 2 points. " + std::to_string(v) + " " + std::to_string(sqrt(G*(dy + sqrt(dx*dx + dy*dy)))));
  }
  auto const A = v2 / G / dx;
  auto const B = sqrt(D) / G / dx;
  auto const tg = A - B;
  return tg;
}

/// Calculate path to the surface of an inclined plane
/// x0, y0 - initial coordinates of the particle
/// vx, vy - initial velocity of the particle
/// xm, ym - coordinates on a surface of a plane
/// theta - angle the plane makes with the horizon
Path calcPath(double x0, double y0, double vx, double vy, double xm, double ym, double theta) {
  auto const c = cos(theta);
  auto const s = sin(theta);

  auto t1 = calc_t1(vx, vy, G, c, s, x0, y0, xm, ym);


  double const x1 = x0 + vx*t1;
  double const y1 = y0 + t1 * (vy - G * t1 / 2.0);

  // incident velocity
  double const v0x = vx;
  double const v0y = vy - G * t1;

  // reflected velocity
  double const v1x = c * (c * v0x + s * v0y) + s * (c * v0y - s * v0x);
  double const v1y = -c * (c * v0y - s * v0x) + s * (c * v0x + s * v0y);

  return Path{t1, x0, y0, vx, vy, x1, y1, v1x, v1y};
}

/// Calculate path to a vertical plane.
/// x0, y0 - initial coordinates of the particle
/// vx, vy - initial velocity of the particle
/// xm - x-coordinate of a vertical plane where path ends.
Path calcPath(double x0, double y0, double vx, double vy, double xm, bool reflect=false) {

  auto t1 = (xm - x0) / vx;

  if (t1 < 0) {
    throw TrajectoryError("Cannot go back in time! " + std::to_string(xm) + " " + std::to_string(x0));
  }

  double const x1 = xm;
  double const y1 = y0 + t1 * (vy - G * t1 / 2.0);

  // incident velocity
  double const v0x = vx;
  double const v0y = vy - G * t1;

  // reflected velocity
  double const v1x = reflect ? - v0x : v0x;
  double const v1y = v0y;

  return Path{t1, x0, y0, vx, vy, x1, y1, v1x, v1y};
}

Path calcPath(const Path &p, double xm, double ym, double theta) {
  return calcPath(p.x1, p.y1, p.v1x, p.v1y, xm, ym, theta);
}

Path calcPath(const Path &p, double xm) {
  return calcPath(p.x1, p.y1, p.v1x, p.v1y, xm);
}

Path calcPathThrough2Points(double x0, double y0, double x1, double y1, double v) {
  auto const tg = calcAngleToPassBetween2Points(x0, y0, x1, y1, v);
  auto const vx = v / sqrt(tg*tg + 1.0);
  auto const vy = vx * tg;
  auto path = calcPath(x0, y0, vx, vy, x1);
  if (fabs(path.y1 - y1) > 1e-4) {
    std::cerr << "Error in calculating path throught 2 points: " + std::to_string(path.y1) << " " << std::to_string(y1) << std::endl;
    std::cerr << "tg: " << tg << ' ' << atan(tg) / M_PI * 180.0 << std::endl;
    std::cerr << "r0: " << x0 << ' ' << y0 << std::endl;
    std::cerr << "v0: " << vx << ' ' << vy << ' ' << sqrt(vx*vx + vy*vy) << ' ' << v << std::endl;
    throw TrajectoryError("Error in calculating path throught 2 points: " + std::to_string(path.y1) + " " + std::to_string(y1));
  }
  return path;
}

//----------------------------------------------------------------------------------------------
enum class DetectorOrientation {Reflection, Transmission};

struct Setup {
  Setup(double _x1, double _y1, double _x2, double _y2, double _ys,
        double _theta, double _xd, double _yd, double _thetaD, bool _hasSample)
      : x1(_x1), y1(_y1), x2(_x2), y2(_y2), ys(_ys), theta(_theta), xd(_xd),
        yd(_yd), thetaD(_thetaD), hasSample(_hasSample) {}
  Setup(Instrument const &instrument, double const slitTheta, double const totalTheta, double const ySample, DetectorOrientation const detOrient) {
    IComponent_const_sptr slit1 = instrument.getComponentByName("slit1");
    IComponent_const_sptr slit2 = instrument.getComponentByName("slit2");
    IComponent_const_sptr linearDetector = instrument.getComponentByName("LinearDetector");
  
    auto tg = tan(slitTheta);
    x1 = slit1->getPos().Z();
    y1 = fabs(x1) * tg;
    x2 = slit2->getPos().Z();
    y2 = fabs(x2) * tg;
    ys = ySample;
    theta = totalTheta - slitTheta;
    auto r = linearDetector->getPos().Z();
    xd = r * cos(slitTheta);
    yd = r * sin(slitTheta);
    if (detOrient == DetectorOrientation::Reflection) {
      thetaD = slitTheta - M_PI/2;
      hasSample = true;
    } else if (detOrient == DetectorOrientation::Transmission) {
      thetaD =  M_PI/2 - slitTheta;
      yd = - yd;
      hasSample = false;
    } else {
      throw std::runtime_error("Unknown detector orientation");
    }
  }
  double x1;
  double y1;
  double x2;
  double y2;
  double ys;
  double theta;
  double xd;
  double yd;
  double thetaD;
  bool hasSample = true;
  double getSlit1Low() const {return ys + x1*tan(theta);}
  double getSlit2Low() const {return ys + x2*tan(theta);}
  double getSlit1LowestShift(double s) const {
    if (y1 + s < getSlit1Low()) {
      return getSlit1Low() * 0.9 - y1;
    } else {
      return s;
    }
  }
  double getSlit2LowestShift(double s) const {
    if (y2 + s < getSlit2Low()) {
      return getSlit2Low() * 0.9 - y2;
    } else {
      return s;
    }
  }
};

std::ostream &operator<<(std::ostream &ostr, Setup const &s) {
  ostr << "Setup:\n";
  ostr << "    slit1: (" << s.x1 << ',' << s.y1 << ")\n";
  ostr << "    slit2: (" << s.x2 << ',' << s.y2 << ")\n";
  ostr << "       ys: " << s.ys << '\n';
  ostr << "    theta: " << s.theta / M_PI * 180 << '\n';
  ostr << "      det: (" << s.xd << ',' << s.yd << ")\n";
  ostr << "   thetaD: " << s.thetaD / M_PI * 180 << '\n';
  return ostr;
}

//----------------------------------------------------------------------------------------------

Trajectory calcTwoSlitTrajectory(Setup const &s, double v, double dy1 = 0.0,
                                 double dy2 = 0.0,
                                 bool can_miss_sample = true) {
  Trajectory traj;
  if (s.y1 + dy1 < s.getSlit1Low()) {
    throw TrajectoryError("Point in slit 1 is below sample surface.");
  }
  if (s.y2 + dy2 < s.getSlit2Low()) {
    throw TrajectoryError("Point in slit 2 is below sample surface.");
  }
  auto calcPath3 = [](Path const &path, double xd, double yd, double angle) {
    if (angle == 0.0) {
      return calcPath(path, xd);
    } else {
      return calcPath(path, xd, yd, angle);
    }
  };
  auto path1 = calcPathThrough2Points(s.x1, s.y1 + dy1, s.x2, s.y2 + dy2, v);
  traj.push_back(path1);
  auto path3 = calcPath3(path1, s.xd, s.yd, s.thetaD);
  if (s.hasSample) {
    auto path2 = calcPath(path1, 0, s.ys, s.theta);
    if (path3.displacementLength2() < path2.displacementLength2()) {
      if (!can_miss_sample) {
        std::string msg("Particle missed the sample: x1(" +
                        std::to_string(path2.x1) + ") > xd(" +
                        std::to_string(s.xd) + ")\n");
        msg.append("Slit heights: " + std::to_string(s.y1 + dy1) + ", " +
                   std::to_string(s.y2 + dy2));
        throw TrajectoryError(msg);
      }
    } else {
      path3 = calcPath3(path2, s.xd, s.yd, s.thetaD);
      traj.push_back(path2);
    }
  }
  traj.push_back(path3);
  return traj;
}

//----------------------------------------------------------------------------------------------
double getTheta(Setup const &s, Trajectory const& traj) {
  auto const& path = traj.back();
  auto const theta = (atan(path.vy / path.vx) - s.theta);
  return theta / M_PI * 180.0;
}

//----------------------------------------------------------------------------------------------
MatrixWorkspace_sptr makeTwoSlitTrajectory(Setup const &s, double v, bool can_miss_sample=true) {
  auto traj = calcTwoSlitTrajectory(s, v, 0, 0, can_miss_sample);
  std::cerr << "y3 " << traj.back().y1 << std::endl;
  std::cerr << "th " << getTheta(s, traj) << std::endl;

  size_t nSpec = traj.size() + 2;
  size_t const nBins = 100;
  auto out =
      WorkspaceFactory::Instance().create("Workspace2D", nSpec, nBins, nBins);
  size_t i = 0;
  fillSpectrum(*out, i++, -1.0, 1.0, Plane{0, s.ys, s.theta});
  if (s.thetaD != 0) {
    fillSpectrum(*out, i++, -1.0, 1.0, Plane{s.xd, s.yd, s.thetaD});
  } else {
    fillSpectrum(*out, i++, -1.0, 1.0, Vertical{s.xd});
  }
  for(auto &&path: traj) {
    fillSpectrum(*out, i, 0.0, path.t, path);
    ++i;
  }
  return out;
}

struct Reflection {
  double x;
  double y;
  double theta;
};

//----------------------------------------------------------------------------------------------
Reflection calcReflections(Setup const &s, double v, double dy0=0.0, double dy1=0.0, bool can_miss_sample=false) {
  auto traj = calcTwoSlitTrajectory(s, v, dy0, dy1, can_miss_sample);
  auto const& path = traj.back();
  auto const trueTheta = (atan(path.vy / path.vx) - s.theta);
  return Reflection{path.x1, path.y1, trueTheta};
}

//----------------------------------------------------------------------------------------------
double getDetectorTangent(Setup const &s, double yd) {
  auto detY = yd * fabs(sin(s.thetaD)) + s.yd;
  auto detX = yd * fabs(cos(s.thetaD)) + s.xd;
  return detY/detX;
}

//----------------------------------------------------------------------------------------------
double calcY(Setup const &s, Reflection const &refl, bool isCorrected=false) {
  auto sinThetaD = fabs(sin(s.thetaD));
  auto ry = refl.y;
  auto rx = refl.x;
  if (isCorrected) {
    double const tgd = tan(s.thetaD);
    double const vrr = tan(refl.theta + s.theta);
    rx = (tgd*s.xd - s.yd)/(tgd - vrr);
    ry = vrr * rx;
  }
  double y = sinThetaD > 0.0 ? (ry - s.yd) / sinThetaD : rx - s.xd;
  return y;
}

//----------------------------------------------------------------------------------------------
std::pair<double, double> calcYTheta(Setup const &s, Reflection const &refl) {
  return std::make_pair(calcY(s, refl), refl.theta);
}

//----------------------------------------------------------------------------------------------
void makeDebugWorkspaces(Setup const &s, MatrixWorkspace_sptr inputWS, MatrixWorkspace_sptr &xWS, MatrixWorkspace_sptr &yWS) {
  auto n_spec = inputWS->getNumberHistograms();
  auto n_bins = inputWS->blocksize();
  xWS = inputWS->clone();
  yWS = inputWS->clone();
  double const detLength = .05;
  auto xAxis1 = new NumericAxis(n_spec);
  auto yAxis1 = new NumericAxis(n_spec);
  auto dx = fabs(detLength * cos(s.thetaD) / (n_spec - 1));
  auto x0 = s.xd - fabs(detLength * cos(s.thetaD) / 2);
  auto dy = fabs(detLength * sin(s.thetaD) / (n_spec - 1));
  auto y0 = s.yd - fabs(detLength * sin(s.thetaD) / 2);
  for(size_t i = 0; i < n_spec; ++i) {
    xAxis1->setValue(i, x0 + double(i) * dx);
    yAxis1->setValue(i, y0 + double(i) * dy);
    xWS->mutableY(i).assign(n_bins, 0.0);
    yWS->mutableY(i).assign(n_bins, 0.0);
  }
  xWS->replaceAxis(1, xAxis1);
  yWS->replaceAxis(1, yAxis1);
}

//----------------------------------------------------------------------------------------------
struct SimulationParams {
  double ys;
  double theta;
  double slit1Gap;
  double slit2Gap;
  DetectorOrientation detOrient = DetectorOrientation::Reflection;
  double sigma = 0.001;
  double rebinStart = 1.0;
  double rebinEnd = 16.0;
  double rebinStep = 0.1;
  bool isCorrected = false;
};

MatrixWorkspace_sptr makeTwoSlitSimulation(SimulationParams const &params) {
  auto loader = FrameworkManager::Instance().createAlgorithm("LoadEmptyInstrument");
  loader->setChild(true);
  loader->setProperty("Filename", "D:\\Work\\mantid_stuff\\Reflectometry\\FakeFacility\\instruments\\REFL_Definition_0.xml");
  loader->setProperty("OutputWorkspace", "dummy");
  loader->execute();
  MatrixWorkspace_sptr ws = loader->getProperty("OutputWorkspace");
  ws->getAxis(0)->setUnit("Wavelength");
  for(size_t i = 0; i < ws->getNumberHistograms(); ++i) {
    auto &X = ws->mutableX(i);
    X.front() = params.rebinStart;
    X.back() = params.rebinEnd;
    ws->mutableY(i).front() = 0.0;
  }
  auto rebin = FrameworkManager::Instance().createAlgorithm("Rebin");
  rebin->setChild(true);
  rebin->setProperty("InputWorkspace", ws);
  rebin->setProperty("Params", std::vector<double>{params.rebinStart, params.rebinStep, params.rebinEnd});
  rebin->setProperty("OutputWorkspace", "dummy");
  rebin->execute();
  ws = rebin->getProperty("OutputWorkspace");

  Instrument_const_sptr instrument = ws->getInstrument();
  double const slitTheta = 2.3 / 180.0 * M_PI;
  Setup const s(*instrument, slitTheta, params.theta, params.ys, params.detOrient);

  size_t const ny = ws->getNumberHistograms();
  std::vector<double> yd(ny);
  auto detInfo = ws->detectorInfo();
  double yMin = std::numeric_limits<double>::max();
  double yMax = -yMin;
  for(size_t i = 0; i < ny; ++i) {
    if (!detInfo.isMonitor(i)) {
      auto y = detInfo.position(i).Y();
      yd[i] = y;
      if (y < yMin) yMin = y;
      if (y > yMax) yMax = y;
    }
  }

  double yc, thetaC;
  std::tie(yc, thetaC) = calcYTheta(s, calcReflections(s, lam_to_speed, 0, 0, true));
  auto shift = yc - (yMin + yMax) / 2.0;
  std::transform(yd.begin(), yd.end(), yd.begin(), [shift](double a){return a + shift;});
  std::cerr << "yc: " << yc << " theta: " << thetaC / M_PI * 180.0 << std::endl;
  std::cerr << s << std::endl;

  auto wavelengths = ws->points(0);
  for(size_t j = 0; j < wavelengths.size(); ++j) {
    auto lam = wavelengths[j];
    double const v = lam_to_speed / lam;
    for (size_t k = 0; k < 10000; ++k) {
      auto const a = std::uniform_real_distribution<double>(-params.slit1Gap / 2, params.slit1Gap / 2)(rng);
      double b;
      if (params.sigma > 0) {
        auto const x = Kernel::normal_distribution<double>(0.0, std::abs(params.sigma))(rng);
        b = a - x;
        if (b < -params.slit2Gap / 2 || b > params.slit2Gap / 2) continue;
      } else {
        b = std::uniform_real_distribution<double>(-params.slit2Gap / 2, params.slit2Gap / 2)(rng);
      }
      double y;
      try {
        y = calcY(s, calcReflections(s, v, a, b, true), params.isCorrected);
      } catch(TrajectoryError&) {
        continue;
      }
      double yMin = std::numeric_limits<double>::max();
      size_t ii = ny;
      for(size_t i = 0; i < ny; ++i) {
        auto dy = fabs(yd[i] - y);
        if (dy < yMin) {
          yMin = dy;
          ii = i;
        }
      }
      if (ii < ny) {
        ws->mutableY(ii)[j] += 1.0;
      }
    }
  }
  return ws;
}

//----------------------------------------------------------------------------------------------
Setup getSlitPos(Setup const &s, double v, double vrr) {
  auto K = G / (2.0 * v * v);
  auto tg = tan(s.theta);
  auto const vr = (tg*tg*vrr + 2.0*tg - vrr) / (-tg*tg + 2.0*tg*vrr + 1.0);
  auto const tg2 = pow(vr - tg, 2);
  auto A = -K * s.x1 * s.x1 + 2 * K * s.x1 * s.x2 - K * s.x2 * s.x2 + s.y1 - tg * (s.x1 - s.x2);
  auto D = tg2 - 4 * K * s.y1 + 4 * K * tg * s.x1;
  if (D < 0.0) {
    throw TrajectoryError("Cannot get position in slit 2");
  }
  auto B = sqrt(D) * (s.x1 - s.x2);
  auto ds2_1 = A - B;
  auto ds2_2 = A + B;
  auto ds2 = fabs(ds2_1 - s.y2) < fabs(ds2_2 - s.y2) ? ds2_1 : ds2_2;
  Setup s1(s);
  if (ds2 > s.getSlit2Low()) {
    s1.y2 = ds2;
  } else {
    auto a = (K*s.x1*s.x1 - K*s.x2*s.x2 + s.y1 - ds2)/(s.x1 - s.x2);
    auto y0 = (s.x1*(-K*s.x1*s.x2 + K*s.x2*s.x2 + ds2) - s.x2*s.y1)/(s.x1 - s.x2);
    auto xr = ((a - tg + sqrt(a*a - 2*a*tg + 4*K*y0 + tg*tg))/(2*K) + s.x1) / 2;
    s1.x2 = xr;
    s1.y2 = y0 + xr*(a - K*xr);
  }
  return s1;
}

MatrixWorkspace_sptr applyGravityCorrection(MatrixWorkspace_sptr inputWS, double theta, double ys, std::string const &detctorName, DetectorOrientation detOrient) {
  if (inputWS->getAxis(0)->unit()->unitID() != "Wavelength") {
    throw std::runtime_error(
        "Input workspace is expected to be in wavelength, found " +
        inputWS->getAxis(0)->unit()->unitID());
  }

  auto ws = MatrixWorkspace_sptr(inputWS->clone());
  Instrument_const_sptr instrument = ws->getInstrument();
  double const slitTheta = 2.3 / 180.0 * M_PI;
  Setup const s(*instrument, slitTheta, theta, ys, detOrient);
  auto r = sqrt(s.xd*s.xd + s.yd*s.yd);
  MatrixWorkspace_sptr debugXWS;
  MatrixWorkspace_sptr debugYWS;
  makeDebugWorkspaces(s, inputWS, debugXWS, debugYWS);
  AnalysisDataService::Instance().addOrReplace("debugXcorr", debugXWS);
  AnalysisDataService::Instance().addOrReplace("debugYcorr", debugYWS);

  auto const detInfo = ws->detectorInfo();
  size_t nd = 0;
  double yMin = std::numeric_limits<double>::max();
  double yMax = - yMin;
  size_t n_spec = ws->getNumberHistograms();
  std::vector<double> yd(n_spec);
  std::vector<bool> isGood(n_spec, false);
  for(size_t index = 0; index < n_spec; ++index) {
    if (detInfo.isMonitor(index)) continue;
    auto const pos = detInfo.position(index);
    if (fabs(pos.Z() - r) > 1e-4) continue;
    auto const y = pos.Y();
    yd[index] = y;
    isGood[index] = true;
    if (y < yMin) yMin = y;
    if (y > yMax) yMax = y;
    ++nd;
  }
  for (size_t i = 0; i < n_spec; ++i) {
    if (!isGood[i])
      continue;
    auto &Y = ws->mutableY(i);
    Y.assign(Y.size(), 0.0);
  }
  auto dy = (yMax - yMin) / (nd - 1);

  std::cerr << s << std::endl;
  double yc, thetaC;
  std::tie(yc, thetaC) = calcYTheta(s, calcReflections(s, lam_to_speed));
  std::cerr << "Central det: " << yc << ' '  << (yc - yMin) / dy << ' ' << thetaC / M_PI * 180.0 << std::endl;
  {
    auto const shift = yc - (yMin + yMax) / 2.0;
    auto alg =
        FrameworkManager::Instance().createAlgorithm("MoveInstrumentComponent");
    alg->setChild(true);
    alg->setProperty("Workspace", ws);
    alg->setProperty("ComponentName", "LinearDetector");
    alg->setProperty("Y", shift);
    alg->execute();
    yMin += shift;
    yMax += shift;
    std::transform(yd.begin(), yd.end(), yd.begin(), [shift](double a){return a + shift;});
  }
  auto const shift = - dy / 2;
  std::transform(yd.begin(), yd.end(), yd.begin(), [shift](double a){return a + shift;});
  yd.push_back(yd.back() + dy);
  std::cerr << yMin << ' ' << yMax << std::endl;
  std::cerr << yd.front() << ' ' << yd.back() << std::endl;

  bool debug = false;

  for(size_t i = 0; i < n_spec; ++i) {
    if (!isGood[i]) continue;
    if (i == n_spec - 1) debug = true;
    auto points = ws->points(i);
    auto &Y = ws->mutableY(i);
    auto y0_lo = yd[i];
    auto y0_hi = yd[i + 1];
    auto tg1_lo = getDetectorTangent(s, y0_lo);
    auto tg1_hi = getDetectorTangent(s, y0_hi);
    if (debug) {
      std::cerr << "y0_lo = " << y0_lo << ' ' << tg1_lo << std::endl;
    }
    for(size_t j = 0; j < points.size(); ++j) {
      auto lam = points[j];
      auto v = lam_to_speed / lam;
      double y_lo, y_hi;
      Setup s1(s);
      try {
        s1 = getSlitPos(s, v, tg1_lo);
        auto refl_lo = calcReflections(s1, v);
        y_lo = calcY(s1, refl_lo);
        if (debug && j == 128) {
          std::cerr << "lam = " << lam << std::endl;
          std::cerr << "y_lo = " << y_lo << std::endl;
        }
        s1 = getSlitPos(s, v, tg1_hi);
        auto refl_hi = calcReflections(s1, v);
        y_hi = calcY(s1, refl_hi);
        size_t k = n_spec;
        for(size_t m = 0; m < n_spec; ++m) {
          auto yd_lo = yd[m];
          auto yd_hi = yd[m + 1];
          double f = 0.0;
          if (y_lo <= yd_lo && yd_lo < y_hi) {
            if (yd_hi < y_hi) {
              f = 1.0;
            } else {
              f = (y_hi - yd_lo) / (yd_hi - yd_lo);
            }
          } else if (y_lo <= yd_hi && yd_hi < y_hi) {
            f = (yd_hi - y_lo) / (yd_hi - yd_lo);
          }
          if (f != 0.0) {
            auto value = inputWS->y(m)[j] * f;
            Y[j] += value;
            try {
              auto xValue = (refl_lo.x + refl_hi.x) / 2;
              //std::cerr << i << ' ' << xValue << ' ' << debugXWS->getAxis(1)->getMin() << " - " << debugXWS->getAxis(1)->getMax() << std::endl;
              size_t ix = debugXWS->getAxis(1)->indexOfValue(xValue);
              debugXWS->mutableY(ix)[j] += value;
              auto yValue = (refl_lo.y + refl_hi.y) / 2;
              size_t iy = debugYWS->getAxis(1)->indexOfValue(yValue);
              debugYWS->mutableY(iy)[j] += value;
            } catch (std::out_of_range &) {
            }
          }
        }
      } catch (TrajectoryError&) {
        continue;
      }
    }
  }
  return ws;
}

} // namespace
//----------------------------------------------------------------------------------------------
void Gravity::init() {
  declareProperty("X0", -2.236, "Source x position in m");
  declareProperty("Y0", 0.0898, "Source y position in m");

  declareProperty("X1", -0.313, "Slit2 x position in m");
  declareProperty("Y1", 0.0126, "Slit2 y position in m");

  declareProperty("Xm1", -10.0, "Mirror x position in m");
  declareProperty("Ym1", 0.0, "Mirror y position in m");
  declareProperty("Theta1", -0.01, "Angle of the mirror in degrees");

  declareProperty("Xd", 3.6, "Position of a linear detector in m");
  declareProperty("Yd", 0.0, "Position of a linear detector in m");
  declareProperty("DetectorName", "LinearDetector", "Name of the detector");
  declareProperty("Theta_D", 0.0, "Position of a linear detector in m");
  declareProperty("Lambda", 0.0, "Particle's wavelength in Angstrom");

  declareProperty("Ys", 0.0, "Vertical position of the sample in m");
  declareProperty("Theta", 0.0, "Angle of the sample in degrees (to the hirizon).");
  declareProperty("Transmission", false, "");

  declareProperty("YSIM0", 0.001, "");
  declareProperty("YSIM1", 0.001, "");
  declareProperty("SIGMA", 0.001, "");
  declareProperty("Corrected", false, "");

  declareProperty(make_unique<WorkspaceProperty<API::MatrixWorkspace>>(
                      "InputWorkspace", "", Direction::Input, PropertyMode::Optional),
                  "The input workspace.");

  declareProperty(make_unique<WorkspaceProperty<API::MatrixWorkspace>>(
                      "OutputWorkspace", "", Direction::Output, PropertyMode::Optional),
                  "The output workspace.");

  declareProperty(make_unique<WorkspaceProperty<API::MatrixWorkspace>>(
                      "PathWorkspace", "", Direction::Output, PropertyMode::Optional),
                  "The output workspace.");
  
  declareProperty(make_unique<WorkspaceProperty<API::MatrixWorkspace>>(
                      "SimulationWorkspace", "", Direction::Output, PropertyMode::Optional),
                  "The output workspace.");
}

//----------------------------------------------------------------------------------------------
void Gravity::exec() {

  double const x0 = getProperty("X0");
  double const y0 = getProperty("Y0");
  double const x1 = getProperty("X1");
  double const y1 = getProperty("Y1");
  double const lam = getProperty("Lambda");
  double const v = lam_to_speed / lam;

  double const xm1 = getProperty("Xm1");
  double const ym1 = getProperty("Ym1");
  double const theta1 = double(getProperty("Theta1")) * M_PI / 180.0;

  double const xd = getProperty("Xd");
  double const yd = getProperty("Yd");
  double const thetaD = double(getProperty("THETA_D")) * M_PI / 180.0;

  double const ys = getProperty("Ys");
  double const theta = double(getProperty("Theta")) * M_PI / 180.0;
  bool isTransmission = getProperty("Transmission");
  auto detOrient = isTransmission ? DetectorOrientation::Transmission : DetectorOrientation::Reflection;

  Setup const s(x0, y0, x1, y1, ys, theta, xd, yd, thetaD, !isTransmission);

  if (!isDefault("InputWorkspace")) {
    MatrixWorkspace_sptr inputWS = getProperty("InputWorkspace");
    std::string detectorName = getProperty("DetectorName");
    MatrixWorkspace_sptr ws = applyGravityCorrection(inputWS, theta, ys, detectorName, detOrient);
    setProperty("OutputWorkspace", ws);
  }

  if (!isDefault("PathWorkspace")) {
    std::cerr << s << std::endl;
    auto pathWS = makeTwoSlitTrajectory(s, v);
    setProperty("PathWorkspace", pathWS);
  }

  if (!isDefault("SimulationWorkspace")) {
    double const ySim0 = getProperty("YSIM0");
    double const ySim1 = getProperty("YSIM1");
    double const sigma = getProperty("SIGMA");
    bool const isCorrected = getProperty("Corrected");
    SimulationParams params;
    params.ys = ys;
    params.theta = theta;
    params.slit1Gap = ySim0;
    params.slit2Gap = ySim1;
    params.detOrient = detOrient;
    params.sigma = sigma;
    params.isCorrected = isCorrected;
    if (!isDefault("Lambda")) {
      params.rebinStart = lam - 0.01;
      params.rebinStep = 0.02;
      params.rebinEnd = lam + 0.01;
    }
    auto ws = makeTwoSlitSimulation(params);
    setProperty("SimulationWorkspace", ws);
  }
}

} // namespace Algorithms
} // namespace Mantid
