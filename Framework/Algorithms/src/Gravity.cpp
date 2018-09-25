#include "MantidAlgorithms/Gravity.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/TextAxis.h"

#include <iostream>

using namespace Mantid::Kernel;
using namespace Mantid::API;

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
double const lam_to_speed = h_plank / m_neutron;
double const G = 9.80665;

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
  double x(double t) const { return t; }
  double y(double t) const { return (t - xm) * tan(theta) + ym; }
};

struct Vertical {
  double xm;
  double x(double t) const { return xm; }
  double y(double t) const { return t; }
};

struct FreeFall {
  double x0, y0, vx, vy;
  double x(double t) const { return x0 + vx * t; }
  double y(double t) const { return y0 + t * (vy - G * t / 2.0); }
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
    throw std::runtime_error("Particle cannot hit the plane.");
  }
  double const r = c*vy - s*vx;
  auto const t1 = (r - sqrt(d))/(c*g);
  auto const t2 = (r + sqrt(d))/(c*g);
  if (t1 < t2 && t1 > 0) {
    return t1;
  } else {
    return t2;
  }
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
    throw std::runtime_error("Particle cannot pass through 2 points.");
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
    throw std::runtime_error("Cannot go back in time! " + std::to_string(xm) + " " + std::to_string(x0));
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
  auto const vx = v / (tg*tg + 1.0);
  auto const vy = vx * tg;
  auto path = calcPath(x0, y0, vx, vy, x1);
  if (fabs(path.y1 - y1) > 1e-5) {
    throw std::runtime_error("Error in calculating path throught 2 points: " + std::to_string(path.y1) + " " + std::to_string(y1));
  }
  return path;
}

std::pair<double, double> calcVs(double lam, double beamDiv) {
  double const v = lam_to_speed / lam * 1e10;
  auto const vx = v * cos(beamDiv);
  auto const vy = v * sin(beamDiv);
  return std::make_pair(vx, vy);
}

std::pair<double, double> calcReflections(double lam, double beamDiv, double x0,
                                          double y0, double xm, double ym,
                                          double theta1, double xd,
                                          double theta) {
  double vx, vy;
  std::tie(vx, vy) = calcVs(lam, beamDiv);
  std::vector<Path> paths;
  Path p1{0, x0, y0, vx, vy, x0, y0, vx, vy};
  try {
    p1 = calcPath(x0, y0, vx, vy, xm, ym, theta1);
  } catch (std::runtime_error &) {
  }
  Path p2;
  try {
    p2 = calcPath(p1, 0, 0, theta);
  } catch (std::runtime_error &) {
    p2 = p1;
  }
  Path p3 = calcPath(p2, xd);
  auto const trueTheta = (atan(p3.vy / p3.vx) - theta) * 180.0 / M_PI;
  return std::make_pair(p3.y1, trueTheta);
}

//----------------------------------------------------------------------------------------------

MatrixWorkspace_sptr makeTwoSlitSimulation(double v, double x0, double y0,
                                           double x1, double y1, double ys,
                                           double theta, double xd, Logger &logger) {
  auto path1 = calcPathThrough2Points(x0, y0, x1, y1, v);
  logger.warning() << "p1 " << path1.v1x << ' ' << path1.v1y << std::endl;
  auto path2 = calcPath(path1, 0, ys, theta);
  logger.warning() << "p2 " << path2.v1x << ' ' << path2.v1y << ' ' << path2.x1 << std::endl;
  Path path3;
  if (path2.x1 > xd) {
    path3 = calcPath(path1, xd);
  } else {
    path3 = calcPath(path2, xd);
  }
  logger.warning() << "p3 " << path3.vx << ' ' << path3.vy << std::endl;
  logger.warning() << "y3 " << path3.y1 << std::endl;

  size_t nSpec = 3;
  size_t const nBins = 100;
  auto out =
      WorkspaceFactory::Instance().create("Workspace2D", nSpec, nBins, nBins);
  fillSpectrum(*out, 0, 0.0, path1.t, path1);
  fillSpectrum(*out, 1, 0.0, path2.t, path2);
  fillSpectrum(*out, 2, 0.0, path3.t, path3);
  return out;
}
//----------------------------------------------------------------------------------------------

MatrixWorkspace_sptr makeTrueTheta(double beamDiv, double x0,
                                   double y0, double xm, double ym,
                                   double theta1, double xd, double theta,
                                   Logger &logger) {
  size_t const nBins = 100;
  auto out =
      WorkspaceFactory::Instance().create("Workspace2D", 3, nBins, nBins);
  auto &X0 = out->mutableX(0);
  auto &Y0 = out->mutableY(0);
  auto &X1 = out->mutableX(1);
  auto &Y1 = out->mutableY(1);
  auto &X2 = out->mutableX(2);
  auto &Y2 = out->mutableY(2);
  auto dt = (13.0 - 1.0) / (nBins - 1);
  double yd, trueTheta;
  for(size_t i = 0; i < nBins; ++i) {
    auto lam = 1.0 + double(i) * dt;
    std::tie(yd, trueTheta) =
        calcReflections(lam, beamDiv, x0, y0, xm, ym, theta1, xd, theta);
    auto j = nBins - 1 - i;
    X0[j] = yd;
    Y0[j] = trueTheta;
    X1[i] = lam;
    Y1[i] = trueTheta;
    X2[i] = lam * (xd - x0) / lam_to_speed * 1e-10 * 1e6;
    Y2[i] = yd;
  }
  return out;
}

//----------------------------------------------------------------------------------------------

MatrixWorkspace_sptr makeSimulation(double lam, double beamDiv, double x0,
                                    double y0, double xm, double ym,
                                    double theta1, double xd, double theta,
                                    Logger &logger) {
  size_t nSpec = 3;
  size_t const nBins = 100;
  Path p1, p2, p3;

  double vx, vy;
  std::tie(vx, vy) = calcVs(lam, beamDiv);

  try {
    p1 = calcPath(x0, y0, vx, vy, xm, ym, theta1);
    ++nSpec;
    logger.warning() << "t1 = " << p1.t << std::endl;
    logger.warning() << "x1 = " << p1.x1 << std::endl;
    logger.warning() << "y1 = " << p1.y1 << std::endl;

    p2 = calcPath(p1, 0, 0, theta);
    ++nSpec;
    logger.warning() << "t2 = " << p2.t << std::endl;
    logger.warning() << "x2 = " << p2.x1 << std::endl;
    logger.warning() << "y2 = " << p2.y1 << std::endl;
  
    p3 = calcPath(p2, xd);
    ++nSpec;
    logger.warning() << "t3 = " << p3.t << std::endl;
    logger.warning() << "x3 = " << p3.x1 << std::endl;
    logger.warning() << "y3 = " << p3.y1 << std::endl;
    auto const trueTheta = (atan(p3.vy / p3.vx) - theta) * 180.0 / M_PI;
    logger.warning() << "th = " << trueTheta << std::endl;
  } catch (std::runtime_error &e) {
    logger.warning() << "Path" << (nSpec - 2) << ": " << e.what() << std::endl;
  }

  auto out =
      WorkspaceFactory::Instance().create("Workspace2D", nSpec, nBins, nBins);
  auto axis = new TextAxis(nSpec);
  if (nSpec > 0) {
    axis->setLabel(0, "Mirror");
    double x1 = nSpec > 3 ? p1.x1 * 2 - x0 : xd;
    fillSpectrum(*out, 0, x0, x1, Plane{xm, ym, theta1});
  }
  if (nSpec > 1) {
    axis->setLabel(1, "Sample");
    double x1 = nSpec > 4 ? p2.x1 * 2 : xd;
    fillSpectrum(*out, 1, 0, x1, Plane{0, 0, theta});
  }
  if (nSpec > 2) {
    axis->setLabel(2, "Detector");
    fillSpectrum(*out, 2, -0.1, 0.5, Vertical{xd});
  }
  if (nSpec > 3) {
    axis->setLabel(3, "Path1");
    fillSpectrum(*out, 3, 0, p1.t, p1);
  }
  if (nSpec > 4) {
    axis->setLabel(4, "Path2");
    fillSpectrum(*out, 4, 0, p2.t, p2);
  }
  if (nSpec > 5) {
    axis->setLabel(5, "Path3");
    fillSpectrum(*out, 5, 0, p3.t, p3);
  }
  out->replaceAxis(1, axis);
  return out;
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
  declareProperty("Lambda", 10.0, "Particle's wavelength in Angstrom");

  declareProperty("Ys", 0.0, "Vertical position of the sample in m");
  declareProperty("Theta", 0.0, "Angle of the sample in degrees (to the hirizon).");
  
  declareProperty(make_unique<WorkspaceProperty<API::MatrixWorkspace>>(
                      "PathWorkspace", "", Direction::Output),
                  "The output workspace.");
  
  //declareProperty(make_unique<WorkspaceProperty<API::MatrixWorkspace>>(
  //                    "ThetaWorkspace", "", Direction::Output),
  //                "The output workspace.");
}

//----------------------------------------------------------------------------------------------
void Gravity::exec() {

  double const x0 = getProperty("X0");
  double const y0 = getProperty("Y0");
  double const x1 = getProperty("X1");
  double const y1 = getProperty("Y1");
  double const lam = getProperty("Lambda");
  double const v = lam_to_speed / lam * 1e10;

  double const xm1 = getProperty("Xm1");
  double const ym1 = getProperty("Ym1");
  double const theta1 = double(getProperty("Theta1")) * M_PI / 180.0;

  double const xd = getProperty("Xd");
  double const ys = getProperty("Ys");
  double const theta = double(getProperty("Theta")) * M_PI / 180.0;

  auto pathWS = makeTwoSlitSimulation(v, x0, y0, x1, y1, ys, theta, xd, g_log);
  setProperty("PathWorkspace", pathWS);


}

} // namespace Algorithms
} // namespace Mantid
