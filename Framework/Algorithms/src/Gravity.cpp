#include "MantidAlgorithms/Gravity.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/TextAxis.h"
#include "MantidAPI/NumericAxis.h"
#include "MantidKernel/MersenneTwister.h"
#include "MantidKernel/PseudoRandomNumberGenerator.h"
#include "MantidKernel/normal_distribution.h"
#include "MantidGeometry/Instrument.h"
#include "MantidGeometry/Instrument/DetectorInfo.h"
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
double const lam_to_speed = h_plank / m_neutron;
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
    throw TrajectoryError("Particle cannot hit the plane.");
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
struct Setup {
  double x0;
  double y0;
  double x1;
  double y1;
  double ys;
  double theta;
  double xd;
  double deltaS;
};

std::ostream &operator<<(std::ostream &ostr, Setup const &s) {
  ostr << "Setup: " << s.x0 << ' ' << s.y0 << ' ' << s.xd ;
  return ostr;
}

//----------------------------------------------------------------------------------------------

Trajectory calcTwoSlitTrajectory(Setup const &s, double v, double dy0=0.0, double dy1=0.0, bool can_miss_sample=true) {
  Trajectory traj;
  if (s.y0 + dy0 - s.x0*tan(s.theta) < 0) throw TrajectoryError("Slit 1 is below 0.");
  if (s.y1 + dy1 - s.x1*tan(s.theta) < 0) throw TrajectoryError("Slit 2 is below 0.");
  auto path1 = calcPathThrough2Points(s.x0, s.y0 + dy0, s.x1, s.y1 + dy1, v);
  traj.push_back(path1);
  auto path2 = calcPath(path1, 0, s.ys, s.theta);
  Path path3;
  if (path2.x1 > s.xd) {
    if (!can_miss_sample) {
      throw TrajectoryError("Particle missed the sample: " + std::to_string(path2.x1) + " > " + std::to_string(s.xd));
    }
    path3 = calcPath(path1, s.xd);
  } else {
    path3 = calcPath(path2, s.xd);
    traj.push_back(path2);
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
MatrixWorkspace_sptr makeTwoSlitTrajectory(Setup const &s, double v, Logger &logger) {
  auto traj = calcTwoSlitTrajectory(s, v);
  logger.warning() << "y3 " << traj.back().y1 << std::endl;
  logger.warning() << "th " << getTheta(s, traj) << std::endl;

  size_t nSpec = traj.size();
  size_t const nBins = 100;
  auto out =
      WorkspaceFactory::Instance().create("Workspace2D", nSpec, nBins, nBins);
  size_t i = 0;
  for(auto &&path: traj) {
    fillSpectrum(*out, i, 0.0, path.t, path);
    ++i;
  }
  return out;
}

//----------------------------------------------------------------------------------------------
std::pair<double, double> calcReflections(Setup const &s, double v, double dy0=0.0, double dy1=0.0) {
  auto traj = calcTwoSlitTrajectory(s, v, dy0, dy1, false);
  auto const& path = traj.back();
  auto const trueTheta = (atan(path.vy / path.vx) - s.theta);
  return std::make_pair(path.y1, trueTheta);
}

//----------------------------------------------------------------------------------------------
double calcLevel(Setup const &setup, double v, double h, double a, double b) {
  auto traj = calcTwoSlitTrajectory(setup, v, a, b, false);
  auto const& path = traj.back();
  return path.y1 - h;
}

//----------------------------------------------------------------------------------------------
std::pair<double, double> calcYTheta(Setup const &s, double v, double ds2=0.0) {

  auto ds20 = - s.deltaS / 2;
  auto const ds21 = s.deltaS / 2;
  if (s.y1 + ds20 < s.x1 * tan(s.theta)) {
    ds20 = 0.9 * (- s.y1 + s.x1 * tan(s.theta));
  }
  double const sigma = 0.001;
  double const da = 1e-5;
  double y0, theta0;
  std::tie(y0, theta0) = calcReflections(s, v, 0, ds2);

  auto const dy1 = (calcLevel(s, v, 0.0, da, ds2) - y0) / da;
  auto const dy2 = (calcLevel(s, v, 0.0, 0, ds2 + da) - y0) / da;
  auto const K = -dy1 / dy2;
  double yd, theta;
  std::tie(yd, theta) = calcReflections(s, v, da, ds2 + K * da);
  double const R = (theta - theta0) / da;
  auto a0 = ds20;
  if (ds2 + K * a0 < ds20) {
    a0 =  (ds20 - ds2) / K;
  }
  auto a1 = ds21;
  if (ds2 + K * a1 > ds21) {
    a1 = (ds21 - ds2) / K;
  }
  auto th0 = theta0 + R * a0;
  auto th1 = theta0 + R * a1;
  // More realistic but doesn't fit by a quadratic
  //auto trueTheta = (th0 + th1) / 2;
  auto trueTheta = theta0;
  return std::make_pair(y0, trueTheta);
}

//----------------------------------------------------------------------------------------------
double findAtLevel(Setup const &setup, double v, double h, double a, double b00) {
  double const tol = 1e-5;
  auto fun = [&](double b){return calcLevel(setup, v, h, a, b);};
  double b0 = a;
  double f0 = fun(b0);
  double b1 = b00;
  double f1 = fun(b1);
  auto b2 = b1;
  while (fabs(f1) > tol) {
    b2 = b1 - f1 *(b1 - b0) / (f1 - f0);
    f0 = f1;
    f1 = fun(b2);
    b0 = b1;
    b1 = b2;
  }
  return b2;
}

MatrixWorkspace_sptr makeYSurface(Setup const &s, double v, Logger &logger) {
  size_t const nBins = 20;
  auto out =
      WorkspaceFactory::Instance().create("Workspace2D", nBins, nBins, nBins);
  auto axis1 = new NumericAxis(nBins);
  out->replaceAxis(1, axis1);
  auto da = (0.1 - -0.1) / (nBins - 1);
  for(size_t i = 0; i < nBins; ++i) {
    auto b = -0.1 + double(i) * da;
    axis1->setValue(i, b);
    auto &X = out->mutableX(i);
    auto &Y = out->mutableY(i);
    for(size_t j = 0; j < nBins; ++j) {
      auto a = -0.1 + double(j) * da;
      X[j] = a;
      Y[j] = calcLevel(s, v, 0.0, a, b);
    }
  }
  return out;
}
//----------------------------------------------------------------------------------------------

MatrixWorkspace_sptr makeTrueTheta(Setup const &s, double lam0, Logger &logger) {
  double const lamStart = lam0 / 10.0;
  double const lamEnd = lam0 * 2;
  size_t const nBins = 20;
  auto out =
      WorkspaceFactory::Instance().create("Workspace2D", 6, nBins, nBins);
  auto axis1 = new TextAxis(nBins);
  axis1->setLabel(0, "Theta/Yd");
  axis1->setLabel(1, "Theta/Lam");
  axis1->setLabel(2, "Yd/Lam");
  axis1->setLabel(3, "b/a");
  axis1->setLabel(4, "Yd/a");
  axis1->setLabel(5, "Theta/a");
  out->replaceAxis(1, axis1);
  auto &X0 = out->mutableX(0);
  auto &Y0 = out->mutableY(0);
  auto &X1 = out->mutableX(1);
  auto &Y1 = out->mutableY(1);
  auto &X2 = out->mutableX(2);
  auto &Y2 = out->mutableY(2);
  auto dt = (lamEnd - lamStart) / (nBins - 1);
  double yd, trueTheta;
  double v0 = 0.0, h0 = 0.0;
  for(size_t i = 0; i < nBins; ++i) {
    auto lam = lamStart + double(i) * dt;
    auto v = lam_to_speed / lam * 1e10;
    std::tie(yd, trueTheta) = calcReflections(s, v);
    trueTheta *= 180.0 / M_PI;
    //auto j = nBins - 1 - i;
    X0[i] = yd;
    Y0[i] = trueTheta;
    X1[i] = lam;
    Y1[i] = trueTheta;
    X2[i] = lam;
    Y2[i] = yd;
    if (lam0 >= lam && (lam0 - lam) <= dt) {
      v0 = v;
      h0 = yd;
    }
  }

  auto &X3 = out->mutableX(3);
  auto &Y3 = out->mutableY(3);
  auto &X4 = out->mutableX(4);
  auto &Y4 = out->mutableY(4);
  auto &X5 = out->mutableX(5);
  auto &Y5 = out->mutableY(5);
  auto da = (0.1 - -0.1) / (nBins - 1);
  double b0 = 0;
  for(size_t i = 0; i < nBins; ++i) {
    auto a = -0.1 + double(i) * da;
    auto b = findAtLevel(s, v0, h0, a, b0);
    b0 = b;
    X3[i] = a;
    Y3[i] = b;
    std::tie(yd, trueTheta) = calcReflections(s, v0, a, b);
    trueTheta *= 180.0 / M_PI;
    X4[i] = a;
    Y4[i] = yd;
    X5[i] = a;
    Y5[i] = trueTheta;
  }
  return out;
}

//----------------------------------------------------------------------------------------------
struct ThetaVsLamY {
  double p0, p1, p2, p3, p4, p5;
  double operator() (double lam, double y) {
    return p0 + p1*y + ((p2 + p3*y) + (p4 + p5*y)*lam)*lam;
  }
};

std::ostream &operator<<(std::ostream &ostr, ThetaVsLamY const &s) {
  ostr << "Func: " << s.p0 << ' ' << s.p1 << ' ' << s.p2 << ' ' << s.p3 << ' ' << s.p4 << ' ' << s.p5 ;
  return ostr;
}

ThetaVsLamY makeThetaFunction(Setup const &setup, double lam1, double lam3) {
  Eigen::Matrix<double, 6, 6> M;
  Eigen::Matrix<double, 6, 1> b;
  double const lam2 = (lam1 + lam3) / 2.0;
  double const s1(0.0), s2(1e-3);
  auto v = [](double lam) {return lam_to_speed / lam * 1e10;};
  double yd, theta;
  size_t row = 0;
  for (auto lam : {lam1, lam2, lam3}) {
    for (auto s : {s1, s2}) {
      std::tie(yd, theta) = calcYTheta(setup, v(lam), s);
      M(row, 0) = 1.0;
      M(row, 1) = yd;
      M(row, 2) = lam;
      M(row, 3) = yd * lam;
      M(row, 4) = lam * lam;
      M(row, 5) = yd * lam * lam;
      b(row) = theta;
      ++row;
    }
  }
  Eigen::Matrix<double, 1, 6> x = M.fullPivLu().solve(b);
  std::cerr << x << std::endl;
  return ThetaVsLamY{x(0), x(1), x(2), x(3), x(4), x(5)};
}

MatrixWorkspace_sptr makeCalcStuff(Setup const &setup, double deltaS, double yd0, double yd1) {
  double const lam1(1.0),  lam3(13.0);
  auto calcTheta = makeThetaFunction(setup, lam1, lam3);

  size_t const ny = 20;
  size_t const nlam = 20;
  auto out =
      WorkspaceFactory::Instance().create("Workspace2D", ny, nlam, nlam);
  auto axis1 = new NumericAxis(ny);
  out->replaceAxis(1, axis1);
  double const dy = (yd1 - yd0) / (ny - 1);
  double const dlam = (lam3 - lam1) / (nlam - 1);
  for(size_t j = 0; j < nlam; ++j) {
    double const lam = 1.0 + dlam * double(j);
    auto &X = out->mutableX(j);
    auto &Y = out->mutableY(j);
    auto &E = out->mutableE(j);
    axis1->setValue(j, lam);
    for(size_t i = 0; i < ny; ++i) {
      X[i] = yd0 + double(i) * dy;
      Y[i] = calcTheta(lam, X[i]) / M_PI * 180.0;
      E[i] = 0;
    }}
  return out;
}

//----------------------------------------------------------------------------------------------
MatrixWorkspace_sptr makeStuff(Setup const &s, double deltaS, double &yd0, double &yd1, bool is_ill = false) {
  yd0 = std::numeric_limits<double>::max();
  yd1 = -yd0;
  size_t const ns = 20;
  size_t const nlam = 20;
  auto out =
      WorkspaceFactory::Instance().create("Workspace2D", ns, nlam, nlam);
  auto axis1 = new NumericAxis(ns);
  out->replaceAxis(1, axis1);

  auto ds20 = - deltaS / 2;
  auto const ds21 = deltaS / 2;
  if (s.y1 + ds20 < s.x1 * tan(s.theta)) {
    ds20 = 0.9 * (- s.y1 + s.x1 * tan(s.theta));
  }
  double const ds = (ds21 - ds20) / (ns - 1);
  double const dlam = (13.0 - 1.0) / (nlam - 1);
  double const sigma = 0.001;
  double const da = ds21 * 1e-3;

  for(size_t j = 0; j < nlam; ++j) {
    double const lam = 1.0 + dlam * double(j);
    auto &X = out->mutableX(j);
    auto &Y = out->mutableY(j);
    auto &E = out->mutableE(j);
    axis1->setValue(j, lam);
    double const v = lam_to_speed / lam * 1e10;
    for(size_t i = 0; i < ns; ++i) {
      double const ds2 = ds20 + ds * double(i);
      try {
        double y0, theta0;
        std::tie(y0, theta0) = calcYTheta(s, v, ds2);
        X[i] = y0;// - s.xd * tan(s.theta);
        Y[i] = theta0 / M_PI * 180.0;
        //E[i] = fabs(th0 - th1) / 2;
      } catch (TrajectoryError &) {
        X[i] = 0.0;
        Y[i] = 0.0;
        E[i] = 0.0;
      }
      if (is_ill) {
        auto k = G / 2 / (v*v);
        auto sy0 = s.y0 - s.x0 * tan(s.theta);
        auto sy1 = s.y1 + ds2 - s.x1 * tan(s.theta);
        auto x0 = (k *(s.x0 * s.x0 - s.x1 * s.x1) + (sy0 - sy1)) / (2*k*(s.x0 - s.x1));
        auto y0 = sy0 + k * pow(s.x0 - x0, 2);
        auto theta_ill = - atan(2*k*x0) / M_PI * 180.0;
        auto dTheta = Y[i] - theta_ill;
        auto theta_ill_0 = - atan((sy0 - sy1)/(s.x0 - s.x1)) / M_PI * 180.0;
        //Y[i] = - atan((sy0 - sy1)/(s.x0 - s.x1)) / M_PI * 180.0; // No correction
        Y[i] = theta_ill; // ILL correction
      }
    }
    if (X.front() > X.back()) {
      for(size_t i = 0; i < ns / 2; ++i) {
        auto k = ns - 1 - i;
        auto tmpX = X[i];
        auto tmpY = Y[i];
        X[i] = X[k];
        Y[i] = Y[k];
        X[k] = tmpX;
        Y[k] = tmpY;
      }
    }
    if (X.front() < yd0) {
      yd0 = X.front();
    }
    if (X.back() > yd1) {
      yd1 = X.back();
    }
  }
  return out;
}

//----------------------------------------------------------------------------------------------
MatrixWorkspace_sptr makeTwoSlitSimulation(Setup const &s, double const yd0, double const yd1) {
  double deltaS = s.deltaS;
  double const sigma = 0.01;
  size_t const ny = 300;
  size_t const nlam = 20;
  auto out =
      WorkspaceFactory::Instance().create("Workspace2D", nlam, ny, ny);
  auto axis1 = new NumericAxis(nlam);
  out->replaceAxis(1, axis1);

  double const dy = (yd1 - yd0) / (ny - 1);

  double const dlam = (13.0 - 1.0) / (nlam - 1);
  for(size_t j = 0; j < nlam; ++j) {
    double const lam = 1.0 + dlam * double(j);
    double const v = lam_to_speed / lam * 1e10;
    auto &X = out->mutableX(j);
    auto &Y = out->mutableY(j);
    auto &E = out->mutableE(j);
    std::vector<double> N(E.size());
    axis1->setValue(j, lam);
    for(size_t i = 0; i < ny; ++i) {
      X[i] = yd0 + double(i) * dy;// - s.xd * tan(s.theta);
      Y[i] = 0;
      E[i] = 0;
    }
    for (size_t k = 0; k < 100000; ++k) {
      //auto const b = std::uniform_real_distribution<double>(-deltaS / 2, deltaS / 2)(rng);
      auto const a = std::uniform_real_distribution<double>(-deltaS / 2, deltaS / 2)(rng);
      auto const x = Kernel::normal_distribution<double>(0.0, std::abs(sigma))(rng);
      auto const b = a - x;
      if (b < -deltaS / 2 || b > deltaS / 2) continue;
      double y, theta;
      try {
        std::tie(y, theta) = calcReflections(s, v, a, b);
      } catch(TrajectoryError&) {
        continue;
      }
      theta /= M_PI / 180.0;
      auto const i = static_cast<size_t>((y - yd0) / dy);
      if (i < ny) {
        N[i] += 1.0;
        Y[i] += theta;
      }
    }
    for(size_t i = 0; i < ny; ++i) {
      if (N[i] > 0)
        Y[i] /= N[i];
      E[i] = atan(X[i]/s.xd)/2 / M_PI * 180.0;
    }
  }
  return out;
}
//----------------------------------------------------------------------------------------------
MatrixWorkspace_sptr applyGravityCorrection(MatrixWorkspace_sptr inputWS, double theta, double deltaS) {
  if (inputWS->getAxis(0)->unit()->unitID() != "Wavelength") {
    throw std::runtime_error(
        "Input workspace is expected to be in wavelength, found " +
        inputWS->getAxis(0)->unit()->unitID());
  }
  auto ws = MatrixWorkspace_sptr(inputWS->clone());
  Instrument_const_sptr instrument = ws->getInstrument();
  IComponent_const_sptr slit1 = instrument->getComponentByName("slit1");
  IComponent_const_sptr slit2 = instrument->getComponentByName("slit2");
  IComponent_const_sptr linearDetector = instrument->getComponentByName("LinearDetector");
  auto tg = tan(theta);
  auto x0 = slit1->getPos().Z();
  auto y0 = fabs(x0) * tg;
  auto x1 = slit2->getPos().Z();
  auto y1 = fabs(x1) * tg;
  double const ys = 0.0;
  auto xd = linearDetector->getPos().Z();
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
    if (pos.Z() != xd) continue;
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
  std::cerr << nd << ' ' << dy << std::endl;
  Setup const s{x0, 0, x1, 0, ys, theta, xd, deltaS};
  std::cerr << s << std::endl;
  auto calcTheta = makeThetaFunction(s, 0.1, 16.0);
  std::cerr << calcTheta << std::endl;
  for(size_t i = 0; i < n_spec; ++i) {
    if (!isGood[i]) continue;
    auto points = ws->points(i);
    //auto &Y = ws->mutableY(i);
    auto const &Y0 = inputWS->y(i);
    auto y0 = yd[i];
    for(size_t j = 0; j < points.size(); ++j) {
      auto lam = points[j];
      auto twoTheta = 2.0*calcTheta(lam, y0);
      double y = xd * tan(twoTheta);
      //std::cerr << yd[i] << ' ' << twoTheta / M_PI * 180.0 / 2 << ' ' << detInfo.twoTheta(i) / M_PI * 180.0 / 2 << std::endl;
      double dyd = fabs(y0 - y);
      size_t k = i + 1;
      while(k < n_spec && dyd > fabs(y0 - yd[k])){
        //std::cerr << k << ' ' << dyd << ' ' << fabs(y0 - yd[k]) << std::endl;
        ++k;
      }
      if (k != n_spec && dyd <= fabs(y0 - yd[k])) {
        ws->mutableY(k)[j] = Y0[j];
        //std::cerr << y0 << ' ' << y << ' ' << dyd << ' ' << yd[k] << std::endl;
      } else if (i > 0) {
        k = i - 1;
        while(k >= 1 && dyd > fabs(y0 - yd[k])) {
          --k;
        }
        if (dyd <= fabs(y0 - yd[k])) {
          ws->mutableY(k)[j] = Y0[j];
        } else {
          ws->mutableY(i)[j] = Y0[j];
        }
      } else {
        ws->mutableY(i)[j] = Y0[j];
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
  declareProperty("DS", 0.01, "Slit2 width in m");

  declareProperty("Xm1", -10.0, "Mirror x position in m");
  declareProperty("Ym1", 0.0, "Mirror y position in m");
  declareProperty("Theta1", -0.01, "Angle of the mirror in degrees");

  declareProperty("Xd", 3.6, "Position of a linear detector in m");
  declareProperty("Lambda", 10.0, "Particle's wavelength in Angstrom");

  declareProperty("Ys", 0.0, "Vertical position of the sample in m");
  declareProperty("Theta", 0.0, "Angle of the sample in degrees (to the hirizon).");

  declareProperty("YSIM0", 0.0, "");
  declareProperty("YSIM1", 0.4, "");

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
                      "ThetaWorkspace", "", Direction::Output, PropertyMode::Optional),
                  "The output workspace.");
  
  declareProperty(make_unique<WorkspaceProperty<API::MatrixWorkspace>>(
                      "YSurfaceWorkspace", "", Direction::Output, PropertyMode::Optional),
                  "The output workspace.");
  
  declareProperty(make_unique<WorkspaceProperty<API::MatrixWorkspace>>(
                      "ILLSurfaceWorkspace", "", Direction::Output, PropertyMode::Optional),
                  "The output workspace.");
  
  declareProperty(make_unique<WorkspaceProperty<API::MatrixWorkspace>>(
                      "SimulationWorkspace", "", Direction::Output, PropertyMode::Optional),
                  "The output workspace.");

  declareProperty(make_unique<WorkspaceProperty<API::MatrixWorkspace>>(
                      "CalcWorkspace", "", Direction::Output, PropertyMode::Optional),
                  "The output workspace.");
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
  double const ds = getProperty("DS");
  double const theta1 = double(getProperty("Theta1")) * M_PI / 180.0;

  double const xd = getProperty("Xd");
  double const ys = getProperty("Ys");
  double const theta = double(getProperty("Theta")) * M_PI / 180.0;

  Setup const s{x0, y0, x1, y1, ys, theta, xd, ds};

  if (!isDefault("InputWorkspace")) {
    MatrixWorkspace_sptr inputWS = getProperty("InputWorkspace");
    auto ws = applyGravityCorrection(inputWS, theta, ds);
    setProperty("OutputWorkspace", ws);
  }

  if (!isDefault("PathWorkspace")) {
    auto pathWS = makeTwoSlitTrajectory(s, v, g_log);
    setProperty("PathWorkspace", pathWS);
  }

  if (!isDefault("ThetaWorkspace")) {
    try {
    auto thetaWS = makeTrueTheta(s, lam, g_log);
    setProperty("ThetaWorkspace", thetaWS);
    } catch (std::runtime_error& e) {
      g_log.warning() << "Cannot calculate ThetaWorkspace: " << e.what() << std::endl;
    }
  }

  if (!isDefault("YSurfaceWorkspace")) {
    double yd0, yd1;
    auto ws = makeStuff(s, ds, yd0, yd1);
    setProperty("YSurfaceWorkspace", ws);
    if (!isDefault("CalcWorkspace")) {
      auto ws = makeCalcStuff(s, ds, yd0, yd1);
      setProperty("CalcWorkspace", ws);
    }
    if (!isDefault("SimulationWorkspace")) {
      double const ySim0 = getProperty("YSIM0");
      double const ySim1 = getProperty("YSIM1");
      auto ws = makeTwoSlitSimulation(s, ySim0, ySim1);
      setProperty("SimulationWorkspace", ws);
    }
    if (!isDefault("ILLSurfaceWorkspace")) {
      auto ws = makeStuff(s, ds, yd0, yd1, true);
      setProperty("ILLSurfaceWorkspace", ws);
    }
  }

}

} // namespace Algorithms
} // namespace Mantid
