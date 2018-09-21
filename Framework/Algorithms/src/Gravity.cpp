#include "MantidAlgorithms/Gravity.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/TextAxis.h"

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

struct Path {
  double t, x0, y0, vx, vy, x1, y1, v1x, v1y;
  double x(double tt) const { return x0 + vx * tt; }
  double y(double tt) const { return y0 + tt * (vy - G * tt / 2.0); }
};

/// Calculate path to the surface of an inclined plane
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
Path calcPath(double x0, double y0, double vx, double vy, double xm) {

  auto t1 = (xm - x0) / vx;

  if (t1 < 0) {
    throw std::runtime_error("Cannot go back in time!");
  }

  double const x1 = xm;
  double const y1 = y0 + t1 * (vy - G * t1 / 2.0);

  // incident velocity
  double const v0x = vx;
  double const v0y = vy - G * t1;

  // reflected velocity
  double const v1x = - v0x;
  double const v1y = v0y;

  return Path{t1, x0, y0, vx, vy, x1, y1, v1x, v1y};
}

Path calcPath(const Path &p, double xm, double ym, double theta) {
  return calcPath(p.x1, p.y1, p.v1x, p.v1y, xm, ym, theta);
}

Path calcPath(const Path &p, double xm) {
  return calcPath(p.x1, p.y1, p.v1x, p.v1y, xm);
}

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

} // namespace
//----------------------------------------------------------------------------------------------
void Gravity::init() {
  declareProperty("X0", -1.0, "Source x position in m");
  declareProperty("Y0", 0.0, "Source y position in m");

  declareProperty("Xm1", 0.0, "Mirror x position in m");
  declareProperty("Ym1", 0.0, "Mirror y position in m");
  declareProperty("Theta1", 0.3, "Angle of the mirror in degrees");

  declareProperty("Xd", 1.0, "Position of a linear detector in m");
  declareProperty("Lambda", 1.0, "Particle's wavelength in Angstrom");
  declareProperty("BeamDiv", 0.00001, "Beam divergence angle in degrees (angle with the horizontal)");

  declareProperty("Theta", 0.3, "Angle of the sample in degrees");
  
  declareProperty(make_unique<WorkspaceProperty<API::MatrixWorkspace>>(
                      "OutputWorkspace", "", Direction::Output),
                  "The output workspace.");
}

//----------------------------------------------------------------------------------------------
void Gravity::exec() {

  double const x0 = getProperty("X0");
  double const y0 = getProperty("Y0");
  double const lam = getProperty("Lambda");
  double const beamDiv = double(getProperty("BeamDiv")) * M_PI / 180.0;

  double const xm1 = getProperty("Xm1");
  double const ym1 = getProperty("Ym1");
  double const theta1 = double(getProperty("Theta1")) * M_PI / 180.0;

  double const xd = getProperty("Xd");
  double const theta = double(getProperty("Theta")) * M_PI / 180.0;

  double const v = lam_to_speed / lam * 1e10;
  auto const vx = v * cos(beamDiv);
  auto const vy = v * sin(beamDiv);

  size_t nSpec = 3;
  size_t const nBins = 100;
  Path p1, p2, p3;

  try {
    p1 = calcPath(x0, y0, vx, vy, xm1, ym1, theta1);
    ++nSpec;
    g_log.warning() << "t1 = " << p1.t << std::endl;
    g_log.warning() << "x1 = " << p1.x1 << std::endl;
    g_log.warning() << "y1 = " << p1.y1 << std::endl;

    p2 = calcPath(p1, 0, 0, theta);
    ++nSpec;
    g_log.warning() << "t2 = " << p2.t << std::endl;
    g_log.warning() << "x2 = " << p2.x1 << std::endl;
    g_log.warning() << "y2 = " << p2.y1 << std::endl;
  
    p3 = calcPath(p2, xd);
    ++nSpec;
    g_log.warning() << "t3 = " << p3.t << std::endl;
    g_log.warning() << "x3 = " << p3.x1 << std::endl;
    g_log.warning() << "y3 = " << p3.y1 << std::endl;
  } catch (std::runtime_error &e) {
    g_log.warning() << "Path" << (nSpec - 2) << ": " << e.what() << std::endl;
  }

  auto out =
      WorkspaceFactory::Instance().create("Workspace2D", nSpec, nBins, nBins);
  auto axis = new TextAxis(nSpec);
  if (nSpec > 0) {
    axis->setLabel(0, "Mirror");
    double x1 = nSpec > 3 ? p1.x1 * 2 - x0 : xd;
    fillSpectrum(*out, 0, x0, x1, Plane{xm1, ym1, theta1});
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

  setProperty("OutputWorkspace", out);


}

} // namespace Algorithms
} // namespace Mantid
