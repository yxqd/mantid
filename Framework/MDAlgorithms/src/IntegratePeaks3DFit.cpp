#include "MantidMDAlgorithms/IntegratePeaks3DFit.h"
#include "MantidDataObjects/MDFramesToSpecialCoordinateSystem.h"
#include "MantidAPI/CommonBinsValidator.h"
#include "MantidAPI/InstrumentValidator.h"
#include "MantidAPI/WorkspaceUnitValidator.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidDataObjects/MDEventWorkspace.h"
#include "MantidDataObjects/MDHistoWorkspace.h"
#include "MantidGeometry/Instrument.h"
#include "MantidKernel/CompositeValidator.h"
#include "MantidGeometry/Crystal/IPeak.h"
#include <boost/math/special_functions/round.hpp>
#include <gsl/gsl_sf_gamma.h> // for factorial
#include <gsl/gsl_linalg.h>   // for SVD
#include <algorithm>
#include <limits>

namespace Mantid {
namespace MDAlgorithms {

using Mantid::Kernel::Direction;
// using Mantid::API::WorkspaceProperty;
using namespace Mantid::DataObjects;
using namespace Mantid::API;
using namespace Mantid::Kernel;
using namespace Mantid::Geometry;

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(IntegratePeaks3DFit)

//----------------------------------------------------------------------------------------------
/**
  * Initialize the algorithm's properties.
  */
void IntegratePeaks3DFit::init() {
  declareProperty(make_unique<WorkspaceProperty<IMDWorkspace> >(
                      "InputWorkspace", "", Direction::Input),
                  "An input Sample MDHistoWorkspace or MDEventWorkspace.");
  declareProperty("Delta", 0.5, "Distance from center to edge of box.");
  declareProperty("GridPoints", 201,
                  "Number of grid points for each dimension of box.");
  declareProperty("NeighborPoints", 10, "Number of points in 5^3 surrounding "
                                        "points above intensity threshold for "
                                        "point to be part of peak.");
  auto fluxValidator = boost::make_shared<CompositeValidator>();
  fluxValidator->add<WorkspaceUnitValidator>("Momentum");
  fluxValidator->add<InstrumentValidator>();
  fluxValidator->add<CommonBinsValidator>();
  auto solidAngleValidator = fluxValidator->clone();

  declareProperty(
      make_unique<WorkspaceProperty<> >("FluxWorkspace", "", Direction::Input,
                                        PropertyMode::Optional, fluxValidator),
      "An optional input workspace containing momentum dependent flux for "
      "normalization.");
  declareProperty(make_unique<WorkspaceProperty<> >(
                      "SolidAngleWorkspace", "", Direction::Input,
                      PropertyMode::Optional, solidAngleValidator),
                  "An optional input workspace containing momentum integrated "
                  "vanadium for normalization "
                  "(a measure of the solid angle).");

  declareProperty(make_unique<WorkspaceProperty<PeaksWorkspace> >(
                      "PeaksWorkspace", "", Direction::Input),
                  "A PeaksWorkspace containing the peaks to integrate.");

  declareProperty(
      make_unique<WorkspaceProperty<PeaksWorkspace> >("OutputWorkspace", "",
                                                      Direction::Output),
      "The output PeaksWorkspace will be a copy of the input PeaksWorkspace "
      "with the peaks' integrated intensities.");
}

//----------------------------------------------------------------------------------------------
/**
 * Execute the algorithm.
 */
void IntegratePeaks3DFit::exec() {
  IMDWorkspace_sptr m_inputWS = getProperty("InputWorkspace");
  Mantid::Kernel::SpecialCoordinateSystem CoordinatesToUse =
      m_inputWS->getSpecialCoordinateSystem();

  /// Peak workspace to integrate
  PeaksWorkspace_sptr inPeakWS = getProperty("PeaksWorkspace");
  const double box = getProperty("Delta");
  const int gridPts = getProperty("GridPoints");
  const int neighborPts = getProperty("NeighborPoints");
  /// Output peaks workspace, create if needed
  PeaksWorkspace_sptr peakWS = getProperty("OutputWorkspace");
  if (peakWS != inPeakWS)
    peakWS = inPeakWS->clone();

  MatrixWorkspace_sptr flux = getProperty("FluxWorkspace");
  MatrixWorkspace_sptr sa = getProperty("SolidAngleWorkspace");

  IMDEventWorkspace_sptr m_eventWS =
      boost::dynamic_pointer_cast<IMDEventWorkspace>(m_inputWS);
  IMDHistoWorkspace_sptr m_histoWS =
      boost::dynamic_pointer_cast<IMDHistoWorkspace>(m_inputWS);
  int npeaks = peakWS->getNumberPeaks();

  auto prog = make_unique<Progress>(this, 0.3, 1.0, npeaks);
  PARALLEL_FOR_IF(Kernel::threadSafe(*peakWS))
  for (int i = 0; i < npeaks; i++) {
    PARALLEL_START_INTERUPT_REGION

    IPeak &p = peakWS->getPeak(i);
    // Get the peak center as a position in the dimensions of the workspace
    V3D pos;
    if (CoordinatesToUse == Mantid::Kernel::QLab) //"Q (lab frame)"
      pos = p.getQLabFrame();
    else if (CoordinatesToUse == Mantid::Kernel::QSample) //"Q (sample frame)"
      pos = p.getQSampleFrame();
    else if (CoordinatesToUse == Mantid::Kernel::HKL) //"HKL"
      pos = p.getHKL();

    MDHistoWorkspace_sptr histoBox;
    if (m_histoWS) {
      histoBox = cropHisto(pos.X(), pos.Y(), pos.Z(), box, m_histoWS);
    } else if (sa && flux) {
      histoBox = normalize(pos.X(), pos.Y(), pos.Z(), box, gridPts, flux, sa,
                           m_eventWS);
    } else {
      histoBox = binEvent(pos.X(), pos.Y(), pos.Z(), box, gridPts, m_eventWS);
    }
    double intensity = 0.0;
    double errorSquared = 0.0;
    integratePeak(neighborPts, histoBox, intensity, errorSquared);
    p.setIntensity(intensity);
    p.setSigmaIntensity(sqrt(errorSquared));
    prog->report();
    PARALLEL_END_INTERUPT_REGION
  }
  PARALLEL_CHECK_INTERUPT_REGION
  // Save the output
  setProperty("OutputWorkspace", peakWS);
}

MDHistoWorkspace_sptr
IntegratePeaks3DFit::normalize(double h, double k, double l, double box,
                               int gridPts, const MatrixWorkspace_sptr &flux,
                               const MatrixWorkspace_sptr &sa,
                               const IMDEventWorkspace_sptr &ws) {
  IAlgorithm_sptr normAlg = createChildAlgorithm("MDNormSCD");
  normAlg->setProperty("InputWorkspace", ws);
  normAlg->setProperty("AlignedDim0",
                       ws->getDimension(0)->getName() + "," +
                           boost::lexical_cast<std::string>(h - box) + "," +
                           boost::lexical_cast<std::string>(h + box) + "," +
                           std::to_string(gridPts));
  normAlg->setProperty("AlignedDim1",
                       ws->getDimension(1)->getName() + "," +
                           boost::lexical_cast<std::string>(k - box) + "," +
                           boost::lexical_cast<std::string>(k + box) + "," +
                           std::to_string(gridPts));
  normAlg->setProperty("AlignedDim2",
                       ws->getDimension(2)->getName() + "," +
                           boost::lexical_cast<std::string>(l - box) + "," +
                           boost::lexical_cast<std::string>(l + box) + "," +
                           std::to_string(gridPts));
  normAlg->setProperty("FluxWorkspace", flux);
  normAlg->setProperty("SolidAngleWorkspace", sa);
  normAlg->setProperty("OutputWorkspace", "mdout");
  normAlg->setProperty("OutputNormalizationWorkspace", "mdnorm");
  normAlg->executeAsChildAlg();
  Workspace_sptr mdout = normAlg->getProperty("OutputWorkspace");
  Workspace_sptr mdnorm = normAlg->getProperty("OutputNormalizationWorkspace");

  IAlgorithm_sptr alg = createChildAlgorithm("DivideMD");
  alg->setProperty("LHSWorkspace", mdout);
  alg->setProperty("RHSWorkspace", mdnorm);
  alg->setPropertyValue("OutputWorkspace", "out");
  alg->execute();
  IMDWorkspace_sptr out = alg->getProperty("OutputWorkspace");
  return boost::dynamic_pointer_cast<MDHistoWorkspace>(out);
}

void IntegratePeaks3DFit::integratePeak(const int neighborPts,
                                        MDHistoWorkspace_sptr out,
                                        double &intensity,
                                        double &errorSquared) {
  // Determine Poision Background Noise, Multivariate Gaussian Distribution   //
  // for Peak and Residual tail signal.                                       //

  // Parameters for Distribution Analysis //

  double eps = 2.2204e-16; // Floating-point relative accuracy
  int do_material = 2;     // Flag to set material files to analyis
                           // do_material = 1 > Silicon
                           // do_material = 2 > Sc
  int n_break_its = 100;   // Max number of iterations for each optimization
  int max_sigma = 30;      // Max STD in Pixels
  int min_sigma = 1;       // Min STD in Pixels
  int max_mu = 10;         // Max deviation from Center in Pixels
  int N_event_thres = 20;  // Minimum number of events to determine signal

  // Begin Program
  int neigh_length_c = 0;
  int neigh_length_m = 0;

  if (do_material == 1) {
    // Silicon file info//
    neigh_length_c =
        3; // Neighborhood size for connection map is 2*neigh_length_c+1
    neigh_length_m =
        3; // Neighborhood size for mean calculation is 2*neigh_length_m+1
  } else if (do_material == 2) {
    // Uncomment for Sc //
    neigh_length_c =
        1; // Neighborhood size for connection map is 2*neigh_length_c+1
    neigh_length_m =
        3; // Neighborhood size for mean calculation is 2*neigh_length_m+1
  } else {
    std::cout << "Material Flag not defined\n";
  }
  std::vector<int> gridPts;
  const size_t dimensionality = out->getNumDims();
  for (size_t i = 0; i < dimensionality; ++i) {
    gridPts.push_back(static_cast<int>(out->getDimension(i)->getNBins()));
  }
  // Record Full signal and background in signal %
  double *num_events = out->getSignalArray();
  double *SqError = out->getErrorSquaredArray();

  // Vectorize Non-zero data and record location //
  std::vector<int> x_vals, y_vals, z_vals, num_events_ref;
  std::vector<double> event_vals, sigma_events;
  double event_tot = 0;
  double event_max = 0;
  int N = gridPts[0];
  double N2 = static_cast<double>(N) / 2.0;
  int N_ind = 0;
  for (int j1 = 0; j1 < gridPts[0]; j1++) {
    for (int j2 = 0; j2 < gridPts[1]; j2++) {
      for (int j3 = 0; j3 < gridPts[2]; j3++) {
        int iPts = j1 + gridPts[0] * (j2 + gridPts[1] * j3);
        // Neither 0, subnormal, infinite, nor NaN
        if (std::isnormal(num_events[iPts])) {
          x_vals.push_back(j1);
          y_vals.push_back(j2);
          z_vals.push_back(j3);
          event_vals.push_back(num_events[iPts]);
          sigma_events.push_back(SqError[iPts]);
          event_tot += num_events[iPts];
          event_max = std::max(event_max, num_events[iPts]);
          num_events_ref.push_back(N_ind);
          N_ind++;
        }
      }
    }
  }

  std::vector<int> event_hist;
  std::vector<int> event_counts;
  int N_non_zeros = 0;
  int number_of_buckets = static_cast<int>(event_max + 0.5);
  for (int j1 = 0; j1 < number_of_buckets; j1++) {
    event_counts.push_back(j1 + 1);
    event_hist.push_back(0);
  }
  for (int j1 = 0; j1 < N_ind; j1++) {
    int bucket = static_cast<int>(event_vals[j1] + 0.5);
    if (bucket > 0)
      event_hist[bucket - 1]++;
  }
  for (int j1 = 0; j1 < number_of_buckets; j1++) {
    if (event_hist[j1] > 0)
      N_non_zeros++;
  }

  double pp_lambda = 0.;

  if (N_ind < N_event_thres) {
    std::cout << "Not enough signal to determine distribution\n";
  } else {
    // Initialize Estimate of Poisson Distribution //
    if (N_non_zeros > 1) {
      pp_lambda = std::max(
          eps, (std::pow((event_hist[1] / gsl_sf_fact(event_counts[1])) /
                             (event_hist[0] / gsl_sf_fact(event_counts[0])),
                         (1.0 / (event_counts[1] - event_counts[0])))));
      std::vector<double> pp_dist;
      for (int j1 = 0; j1 < number_of_buckets; j1++) {
        pp_dist.push_back((std::pow((pp_lambda), event_counts[j1])) *
                          exp(-pp_lambda) / gsl_sf_fact(event_counts[j1]));
      }
      double pp_N = std::inner_product(pp_dist.begin(), pp_dist.end(),
                                       event_hist.begin(), 0.0) /
                    std::inner_product(pp_dist.begin(), pp_dist.end(),
                                       pp_dist.begin(), 0.0);
      double best_pp_err = 0.;
      for (int j1 = 0; j1 < number_of_buckets; j1++) {
        best_pp_err += std::pow((pp_N * pp_dist[j1] - event_hist[j1]), 2);
      }
      best_pp_err /= event_tot;

      double d_pp_lambda = .1 * pp_lambda;

      for (int jpit = 1; jpit < n_break_its; jpit++) {
        double did_mv = 0;
        for (int js = 1; js < 3; js++) {
          double c_pp_lambda = pp_lambda;
          c_pp_lambda += (std::pow((-1), js)) * d_pp_lambda;
          std::vector<double> pp_dist;
          for (int j1 = 0; j1 < number_of_buckets; j1++) {
            pp_dist.push_back((std::pow((c_pp_lambda), event_counts[j1])) *
                              exp(-c_pp_lambda) /
                              gsl_sf_fact(event_counts[j1]));
          }
          double pp_N = std::inner_product(pp_dist.begin(), pp_dist.end(),
                                           event_hist.begin(), 0.0) /
                        std::inner_product(pp_dist.begin(), pp_dist.end(),
                                           pp_dist.begin(), 0.0);
          double c_best_pp_err = 0.;
          for (int j1 = 0; j1 < number_of_buckets; j1++) {
            c_best_pp_err += std::pow((pp_N * pp_dist[j1] - event_hist[j1]), 2);
          }
          c_best_pp_err /= event_tot;
          if (c_best_pp_err < best_pp_err) {
            best_pp_err = c_best_pp_err;
            pp_lambda = c_pp_lambda;
            did_mv = 1;
          }
        }
        if (did_mv == 0)
          d_pp_lambda = d_pp_lambda / 2;
        else
          d_pp_lambda = 1.2 * d_pp_lambda;
        if (d_pp_lambda / pp_lambda < 1e-4)
          break;
      }
    }
  }

  std::vector<double> full_box_mean;
  std::vector<double> full_box_density;
  for (int j0 = 0; j0 < N_ind; j0++) {
    double meanB = 0;
    int num0 = 0;
    int num = 0;
    for (int j1 = x_vals[j0] - neigh_length_m;
         j1 <= x_vals[j0] + neigh_length_m; j1++) {
      if (j1 < 0 || j1 >= gridPts[0])
        continue;
      for (int j2 = y_vals[j0] - neigh_length_m;
           j2 <= y_vals[j0] + neigh_length_m; j2++) {
        if (j2 < 0 || j2 >= gridPts[1])
          continue;
        for (int j3 = z_vals[j0] - neigh_length_m;
             j3 <= z_vals[j0] + neigh_length_m; j3++) {
          if (j3 < 0 || j3 >= gridPts[2])
            continue;
          int iPts = j1 + gridPts[0] * (j2 + gridPts[1] * j3);
          meanB += num_events[iPts];
          num0++;
          if (num_events[iPts] > 0.)
            num++;
        }
      }
    }
    full_box_mean.push_back(meanB / static_cast<double>(num0));
    full_box_density.push_back(
        static_cast<double>(num / static_cast<double>(num0)));
  }
  std::vector<double> conditional_event_vals;
  double conditional_event_vals_tot = 0.0;
  for (int j0 = 0; j0 < N_ind; j0++) {
    if (full_box_mean[j0] >
            pp_lambda + 1.96 * sqrt(pp_lambda /
                                    std::pow((2 * neigh_length_m + 1), 3)) &&
        std::max(fabs(z_vals[j0] - N2),
                 std::max(fabs(x_vals[j0] - N2), fabs(y_vals[j0] - N2))) <=
            max_mu +
                2 * max_sigma) { // Remove points outside of prescribed domain
      conditional_event_vals.push_back(event_vals[j0]);
      conditional_event_vals_tot += event_vals[j0];
    } else
      conditional_event_vals.push_back(0.0);
  }

  // Estimate initial parameters for Gaussian distribution //
  std::vector<double> mu_est(3);

  double muMax = max_mu;
  for (int j0 = 0; j0 < 3; j0++) {
    mu_est[j0] = 0.0;
  }
  for (int j0 = 0; j0 < N_ind; j0++) {
    mu_est[0] +=
        x_vals[j0] * conditional_event_vals[j0] / conditional_event_vals_tot;
    mu_est[1] +=
        y_vals[j0] * conditional_event_vals[j0] / conditional_event_vals_tot;
    mu_est[2] +=
        z_vals[j0] * conditional_event_vals[j0] / conditional_event_vals_tot;
  }
  if (fabs(mu_est[0] - N2 > muMax))
    muMax = fabs(mu_est[0] - N2);
  if (fabs(mu_est[1] - N2 > muMax))
    muMax = fabs(mu_est[1] - N2);
  if (fabs(mu_est[2] - N2 > muMax))
    muMax = fabs(mu_est[2] - N2);
  if (muMax > max_mu) {
    for (int j0 = 0; j0 < 3; j0++) {
      mu_est[j0] = N2;
    }
  }

  DblMatrix covar_mat_est(3, 3);
  for (int jp = 0; jp < N_ind; jp++) {
    covar_mat_est[0][0] += (conditional_event_vals[jp]) *
                           (x_vals[jp] - mu_est[0]) * (x_vals[jp] - mu_est[0]);
    covar_mat_est[1][0] += (conditional_event_vals[jp]) *
                           (x_vals[jp] - mu_est[1]) * (y_vals[jp] - mu_est[0]);
    covar_mat_est[2][0] += (conditional_event_vals[jp]) *
                           (x_vals[jp] - mu_est[2]) * (z_vals[jp] - mu_est[0]);
    covar_mat_est[0][1] += (conditional_event_vals[jp]) *
                           (y_vals[jp] - mu_est[0]) * (x_vals[jp] - mu_est[1]);
    covar_mat_est[1][1] += (conditional_event_vals[jp]) *
                           (y_vals[jp] - mu_est[1]) * (y_vals[jp] - mu_est[1]);
    covar_mat_est[2][1] += (conditional_event_vals[jp]) *
                           (y_vals[jp] - mu_est[2]) * (z_vals[jp] - mu_est[1]);
    covar_mat_est[0][2] += (conditional_event_vals[jp]) *
                           (z_vals[jp] - mu_est[0]) * (x_vals[jp] - mu_est[2]);
    covar_mat_est[1][2] += (conditional_event_vals[jp]) *
                           (z_vals[jp] - mu_est[1]) * (y_vals[jp] - mu_est[2]);
    covar_mat_est[2][2] += (conditional_event_vals[jp]) *
                           (z_vals[jp] - mu_est[2]) * (z_vals[jp] - mu_est[2]);
  }
  for (int j0 = 0; j0 < 3; j0++) {
    for (int j1 = 0; j1 < 3; j1++)
      covar_mat_est[j0][j1] =
          covar_mat_est[j0][j1] / conditional_event_vals_tot;
  }

  gsl_matrix *u_est = gsl_matrix_alloc(3, 3);
  gsl_matrix *v_est = gsl_matrix_alloc(3, 3);
  gsl_vector *d_est = gsl_vector_alloc(3);
  gsl_vector *work = gsl_vector_alloc(3);

  // Need to copy from DblMatrix to gsl matrix

  for (size_t k = 0; k < 3; k++)
    for (size_t l = 0; l < 3; l++)
      gsl_matrix_set(u_est, k, l, covar_mat_est[k][l]);

  gsl_linalg_SV_decomp(u_est, v_est, d_est, work);
  // 2nd column has sign change in gsl from matlab
  double theta_x_est =
      atan2(-gsl_matrix_get(u_est, 2, 1), gsl_matrix_get(u_est, 2, 2));
  double theta_y_est = atan2(-gsl_matrix_get(u_est, 2, 0),
                             sqrt(std::pow(gsl_matrix_get(u_est, 2, 1), 2) +
                                  std::pow(gsl_matrix_get(u_est, 2, 2), 2)));
  double theta_z_est =
      atan2(gsl_matrix_get(u_est, 1, 0), gsl_matrix_get(u_est, 0, 0));

  std::vector<double> sigma(3);
  for (size_t jp = 0; jp < 3; jp++)
    sigma[jp] = gsl_vector_get(d_est, jp);

  for (int jp = 0; jp < 3; jp++) {
    if (sqrt(sigma[jp] / 2.) < min_sigma)
      sigma[jp] = 2.1 * std::pow((min_sigma), 2);
    else if (sqrt(sigma[jp] / 2) > max_sigma)
      sigma[jp] = 1.9 * std::pow((max_sigma), 2);
  }

  std::vector<double> gp_parm;
  for (int jp = 0; jp < 3; jp++) {
    gp_parm.push_back(mu_est[jp]);
  }
  gp_parm.push_back(theta_x_est);
  gp_parm.push_back(theta_y_est);
  gp_parm.push_back(theta_z_est);
  for (int jp = 0; jp < 3; jp++) {
    gp_parm.push_back(sigma[jp]);
  }

  std::vector<double> c_gp_parm = gp_parm;

  DblMatrix X_est(3, 3);
  X_est[0][0] = 1.0;
  X_est[1][0] = 0.0;
  X_est[2][0] = 0.0;
  X_est[0][1] = 0.0;
  X_est[1][1] = std::cos(c_gp_parm[3]);
  X_est[2][1] = -std::sin(c_gp_parm[3]);
  X_est[0][2] = 0.0;
  X_est[1][2] = std::sin(c_gp_parm[3]);
  X_est[2][2] = std::cos(c_gp_parm[3]);

  DblMatrix Y_est(3, 3);
  Y_est[0][0] = std::cos(c_gp_parm[4]);
  Y_est[1][0] = 0.0;
  Y_est[2][0] = std::sin(c_gp_parm[4]);
  Y_est[0][1] = 0.0;
  Y_est[1][1] = 1.0;
  Y_est[2][1] = 0.0;
  Y_est[0][2] = -std::sin(c_gp_parm[4]);
  Y_est[1][2] = 0.0;
  Y_est[2][2] = std::cos(c_gp_parm[4]);

  DblMatrix Z_est(3, 3);
  Z_est[0][0] = std::cos(c_gp_parm[5]);
  Z_est[1][0] = -std::sin(c_gp_parm[5]);
  Z_est[2][0] = 0.0;
  Z_est[0][1] = std::sin(c_gp_parm[5]);
  Z_est[1][1] = std::cos(c_gp_parm[5]);
  Z_est[2][1] = 0.0;
  Z_est[0][2] = 0.0;
  Z_est[1][2] = 0.0;
  Z_est[2][2] = 1.0;

  DblMatrix XYZ_est = Z_est * Y_est * X_est;

  std::vector<double> rot_cord_x;
  std::vector<double> rot_cord_y;
  std::vector<double> rot_cord_z;
  for (int j0 = 0; j0 < N_ind; j0++) {
    std::vector<double> tmp;
    tmp.push_back(x_vals[j0] - c_gp_parm[0]);
    tmp.push_back(y_vals[j0] - c_gp_parm[1]);
    tmp.push_back(z_vals[j0] - c_gp_parm[2]);
    tmp = XYZ_est * tmp;
    rot_cord_x.push_back(tmp[0]);
    rot_cord_y.push_back(tmp[1]);
    rot_cord_z.push_back(tmp[2]);
  }

  std::vector<double> gp_dist;
  for (int j0 = 0; j0 < N_ind; j0++) {
    gp_dist.push_back(exp(-std::pow(rot_cord_x[j0], 2) / c_gp_parm[6]) *
                      exp(-std::pow(rot_cord_y[j0], 2) / c_gp_parm[7]) *
                      exp(-std::pow(rot_cord_z[j0], 2) / c_gp_parm[8]));
  }
  double gp_dist_sum = 0.0;
  for (int j0 = 0; j0 < N_ind; j0++) {
    gp_dist_sum += gp_dist[j0];
  }

  for (int j0 = 0; j0 < N_ind; j0++) {
    gp_dist[j0] /= gp_dist_sum;
  }

  double best_norm = 0;
  double best_norm0 = 0;
  for (int j0 = 0; j0 < N_ind; j0++) {
    best_norm += gp_dist[j0] * conditional_event_vals[j0];
    best_norm0 += gp_dist[j0] * gp_dist[j0];
  }
  best_norm /= best_norm0;

  std::vector<int> gp_events;
  double sum = 0.0;
  for (int j0 = 0; j0 < N_ind; j0++) {
    gp_events.push_back(static_cast<int>(best_norm * gp_dist[j0] + 0.5));
    sum += std::pow((conditional_event_vals[j0] - best_norm * gp_dist[j0]), 2);
  }
  double best_err = sqrt(sum) / event_tot;

  // Estimate Gaussian distribution given Poisson distribution and Extreme Value
  // Distribution //
  std::vector<double> d_gp_parm = gp_parm;
  for (int j0 = 0; j0 < 9; j0++) {
    d_gp_parm[j0] = .01 * fabs(gp_parm[j0]);
    d_gp_parm[j0] =
        std::max(d_gp_parm[j0], .1); // Minimal starting search distance
                                     // of optimization parameters.
  }

  for (int jpit = 0; jpit < n_break_its; jpit++) {
    std::vector<double> did_mv_vec(9);
    for (int jp = 0; jp < 9; jp++)
      did_mv_vec[jp] = 0.0;

    for (int jp = 0; jp < 9; jp++) {
      for (int js = 1; js < 3; js++) {
        c_gp_parm = gp_parm;
        c_gp_parm[jp] = c_gp_parm[jp] + (std::pow((-1), js)) * d_gp_parm[jp];

        int do_opt = 1;
        // Enforce parameter maximum domains %
        if (jp < 3) {
          if (fabs(c_gp_parm[jp] - N2) > max_mu) {
            do_opt = 0;
          }
        }
        if (jp > 5) {
          if (sqrt(c_gp_parm[jp] / 2) < min_sigma ||
              sqrt(c_gp_parm[jp] / 2) > max_sigma) {
            do_opt = 0;
          }
        }

        if (do_opt == 1) {
          X_est[0][0] = 1.0;
          X_est[1][0] = 0.0;
          X_est[2][0] = 0.0;
          X_est[0][1] = 0.0;
          X_est[1][1] = std::cos(c_gp_parm[3]);
          X_est[2][1] = -std::sin(c_gp_parm[3]);
          X_est[0][2] = 0.0;
          X_est[1][2] = std::sin(c_gp_parm[3]);
          X_est[2][2] = std::cos(c_gp_parm[3]);

          Y_est[0][0] = std::cos(c_gp_parm[4]);
          Y_est[1][0] = 0.0;
          Y_est[2][0] = std::sin(c_gp_parm[4]);
          Y_est[0][1] = 0.0;
          Y_est[1][1] = 1.0;
          Y_est[2][1] = 0.0;
          Y_est[0][2] = -std::sin(c_gp_parm[4]);
          Y_est[1][2] = 0.0;
          Y_est[2][2] = std::cos(c_gp_parm[4]);

          Z_est[0][0] = std::cos(c_gp_parm[5]);
          Z_est[1][0] = -std::sin(c_gp_parm[5]);
          Z_est[2][0] = 0.0;
          Z_est[0][1] = std::sin(c_gp_parm[5]);
          Z_est[1][1] = std::cos(c_gp_parm[5]);
          Z_est[2][1] = 0.0;
          Z_est[0][2] = 0.0;
          Z_est[1][2] = 0.0;
          Z_est[2][2] = 1.0;
          XYZ_est = Z_est * Y_est * X_est;

          rot_cord_x.clear();
          rot_cord_y.clear();
          rot_cord_z.clear();
          for (int j0 = 0; j0 < N_ind; j0++) {
            std::vector<double> tmp;
            tmp.push_back(x_vals[j0] - c_gp_parm[0]);
            tmp.push_back(y_vals[j0] - c_gp_parm[1]);
            tmp.push_back(z_vals[j0] - c_gp_parm[2]);
            tmp = XYZ_est * tmp;
            rot_cord_x.push_back(tmp[0]);
            rot_cord_y.push_back(tmp[1]);
            rot_cord_z.push_back(tmp[2]);
          }

          std::vector<double> c_gp_dist;
          for (int j0 = 0; j0 < N_ind; j0++) {
            c_gp_dist.push_back(
                exp(-std::pow(rot_cord_x[j0], 2) / c_gp_parm[6]) *
                exp(-std::pow(rot_cord_y[j0], 2) / c_gp_parm[7]) *
                exp(-std::pow(rot_cord_z[j0], 2) / c_gp_parm[8]));
          }
          double c_gp_dist_sum = 0.0;
          for (int j0 = 0; j0 < N_ind; j0++) {
            c_gp_dist_sum += c_gp_dist[j0];
          }
          // For integration of fit at end
          gp_dist_sum = c_gp_dist_sum;

          for (int j0 = 0; j0 < N_ind; j0++) {
            c_gp_dist[j0] /= c_gp_dist_sum;
          }

          best_norm = 0;
          best_norm0 = 0;
          for (int j0 = 0; j0 < N_ind; j0++) {
            best_norm += c_gp_dist[j0] * conditional_event_vals[j0];
            best_norm0 += c_gp_dist[j0] * c_gp_dist[j0];
          }
          best_norm /= best_norm0;

          std::vector<int> c_gp_events;
          double sum = 0.0;
          for (int j0 = 0; j0 < N_ind; j0++) {
            c_gp_events.push_back(
                static_cast<int>(best_norm * c_gp_dist[j0] + 0.5));
            sum += std::pow(
                (conditional_event_vals[j0] - best_norm * c_gp_dist[j0]), 2);
          }
          double c_best_err = sqrt(sum) / event_tot;
          if (c_best_err < best_err) {
            best_err = c_best_err;
            gp_parm = c_gp_parm;
            gp_dist = c_gp_dist;
            gp_events = c_gp_events;
            did_mv_vec[jp] = 1;
          }
        }
      }
      if (did_mv_vec[jp] == 0) {
        d_gp_parm[jp] = d_gp_parm[jp] / 2.0;
      } else {
        d_gp_parm[jp] = 1.2 * d_gp_parm[jp];
      }
    }

    double max_parm = 0;
    for (int j0 = 0; j0 < 9; j0++) {
      max_parm = std::max(d_gp_parm[j0] / fabs(gp_parm[j0]), max_parm);
    }
    if (max_parm < 1.e-3)
      break;
  }

  // Estimate Poisson distribution given Gaussian distribution //
  std::vector<int> pp_ind, gp_ind;
  for (int j0 = 0; j0 < N_ind; j0++) {
    if (gp_events[j0] == 0) {
      pp_ind.push_back(j0);
    } else {
      gp_ind.push_back(j0);
    }
  }
  event_hist.clear();
  event_counts.clear();
  N_non_zeros = 0;
  number_of_buckets = static_cast<int>(event_max + 0.5);
  for (int j1 = 0; j1 < number_of_buckets; j1++) {
    event_counts.push_back(j1 + 1);
    event_hist.push_back(0);
  }
  for (int j1 = 0; j1 < N_ind; j1++) {
    int bucket = static_cast<int>(event_vals[j1] + 0.5);
    if (bucket > 0)
      event_hist[bucket - 1]++;
  }
  for (int j1 = 0; j1 < number_of_buckets; j1++) {
    if (event_hist[j1] > 0)
      N_non_zeros++;
  }
  pp_lambda = 0.;

  if (N_non_zeros > 1) {
    pp_lambda = std::max(
        eps, (std::pow((event_hist[1] / gsl_sf_fact(event_counts[1])) /
                           (event_hist[0] / gsl_sf_fact(event_counts[0])),
                       (1.0 / (event_counts[1] - event_counts[0])))));
    std::vector<double> pp_dist;
    for (int j1 = 0; j1 < number_of_buckets; j1++) {
      pp_dist.push_back((std::pow((pp_lambda), event_counts[j1])) *
                        exp(-pp_lambda) / gsl_sf_fact(event_counts[j1]));
    }
    double pp_N = std::inner_product(pp_dist.begin(), pp_dist.end(),
                                     event_hist.begin(), 0.0) /
                  std::inner_product(pp_dist.begin(), pp_dist.end(),
                                     pp_dist.begin(), 0.0);
    double best_pp_err = 0.;
    for (int j1 = 0; j1 < number_of_buckets; j1++) {
      best_pp_err += std::pow((pp_N * pp_dist[j1] - event_hist[j1]), 2);
    }
    best_pp_err /= event_tot;

    double d_pp_lambda = .1 * pp_lambda;

    for (int jpit = 1; jpit < n_break_its; jpit++) {
      double did_mv = 0;
      for (int js = 1; js < 3; js++) {
        double c_pp_lambda = pp_lambda;
        c_pp_lambda += (std::pow((-1), js)) * d_pp_lambda;
        std::vector<double> pp_dist;
        for (int j1 = 0; j1 < number_of_buckets; j1++) {
          pp_dist.push_back((std::pow((c_pp_lambda), event_counts[j1])) *
                            exp(-c_pp_lambda) / gsl_sf_fact(event_counts[j1]));
        }
        double pp_N = std::inner_product(pp_dist.begin(), pp_dist.end(),
                                         event_hist.begin(), 0.0) /
                      std::inner_product(pp_dist.begin(), pp_dist.end(),
                                         pp_dist.begin(), 0.0);
        double c_best_pp_err = 0.;
        for (int j1 = 0; j1 < number_of_buckets; j1++) {
          c_best_pp_err += std::pow((pp_N * pp_dist[j1] - event_hist[j1]), 2);
        }
        c_best_pp_err /= event_tot;
        if (c_best_pp_err < best_pp_err) {
          best_pp_err = c_best_pp_err;
          pp_lambda = c_pp_lambda;
          did_mv = 1;
        }
      }
      if (did_mv == 0)
        d_pp_lambda = d_pp_lambda / 2;
      else
        d_pp_lambda = 1.2 * d_pp_lambda;
      if (d_pp_lambda / pp_lambda < 1e-4)
        break;
    }
  }

  std::vector<double> pp_dist;
  for (int j1 = 0; j1 < number_of_buckets; j1++) {
    pp_dist.push_back((std::pow((pp_lambda), event_counts[j1])) *
                      exp(-pp_lambda) / gsl_sf_fact(event_counts[j1]));
  }
  double pp_N =
      std::inner_product(pp_dist.begin(), pp_dist.end(), event_hist.begin(),
                         0.0) /
      std::inner_product(pp_dist.begin(), pp_dist.end(), pp_dist.begin(), 0.0);

  std::vector<int> rp_ind;
  for (int j0 = 0; j0 < N_ind; j0++) {
    if (full_box_mean[j0] >
            pp_lambda + 1.96 * sqrt(pp_lambda /
                                    std::pow((2 * neigh_length_m + 1), 3)) &&
        find(pp_ind.begin(), pp_ind.end(), j0) != pp_ind.end()) {
      rp_ind.push_back(j0);
      std::vector<int>::iterator position =
          std::find(pp_ind.begin(), pp_ind.end(), j0);
      if (position != pp_ind.end())
        pp_ind.erase(position);
    }
  }

  // Estimate Residual Distribution given Poisson distribution and Gaussian
  // distribution //

  /*std::vector<int> connection_map;
  for (int j1 = 0; j1 < N_ind; j1++) {
    if (find(gp_ind.begin(), gp_ind.end(), j1) != gp_ind.end())
      connection_map.push_back(0);
    else
      connection_map.push_back(j1);
  }

  std::vector<double> ind_box;
  for (size_t j0 = 0; j0 < rp_ind.size(); j0++) {
    for (int j1 = x_vals[rp_ind[j0]] - neigh_length_c;
         j1 <= x_vals[rp_ind[j0]] + neigh_length_c; j1++) {
      if (j1 < 0 || j1 >= gridPts[0])
        continue;
      for (int j2 = y_vals[rp_ind[j0]] - neigh_length_c;
           j2 <= y_vals[rp_ind[j0]] + neigh_length_c; j2++) {
        if (j2 < 0 || j2 >= gridPts[1])
          continue;
        for (int j3 = z_vals[rp_ind[j0]] - neigh_length_c;
             j3 <= z_vals[j0] + neigh_length_c; j3++) {
          if (j3 < 0 || j3 >= gridPts[2])
            continue;
          int iPts = j1 + gridPts[0] * (j2 + gridPts[1] * j3);
          ind_box.push_back(num_events_ref[iPts]);
        }
      }
    }
  }
  std::vector<double> cind;
  for (int j0 = 0; j0 < N_ind; j0++) {
    if (full_box_mean[j0] >
        pp_lambda +
            1.96 * sqrt(pp_lambda / std::pow((2 * neigh_length_m + 1), 3))) {
      cind.push_back(static_cast<int>(full_box_mean[j0]));
    }
    for (size_t jp = 1; jp <= rp_ind.size(); jp++) {
      ind_box =
          num_events_ref(max(xyz_vals(rp_ind[jp], 1) - neigh_length_c, 1)
                         : min(xyz_vals(rp_ind[jp], 1) + neigh_length_c, N),
                           max(xyz_vals(rp_ind[jp], 2) - neigh_length_c, 1)
                         : min(xyz_vals(rp_ind[jp], 2) + neigh_length_c, N),
                           max(xyz_vals(rp_ind[jp], 3) - neigh_length_c, 1)
                         : min(xyz_vals(rp_ind[jp], 3) + neigh_length_c, N));
      cind = ind_box(max(xyz_vals(rp_ind[jp], 3) - neigh_length_c, 1)
                     : min(xyz_vals(rp_ind[jp], 3) + neigh_length_c, N));
      cind = ind_box(ind(ind_box ~ = 0));
      cind = cind(find(
          full_box_mean(cind) >
          pp_lambda +
              1.96 * sqrt(pp_lambda / std::pow((2 * neigh_length_m + 1), 3))));

      connection_map(cind( :)) = min(connection_map(cind( :)));
      end

          // Set of points only in the Tail of distribution.  //
          rp_ind = setdiff(find(connection_map == 0), gp_ind);

      // Add points at 5// level of pd //
      neig_ind = [];
      connection_map = zeros(N_ind, 1);
  for
    jp = 1 : length(rp_ind);
  ind_box =
      num_events_ref(max(xyz_vals(rp_ind[jp], 1) - neigh_length_c, 1)
                     : min(xyz_vals(rp_ind[jp], 1) + neigh_length_c, N),
                       ... max(xyz_vals(rp_ind[jp], 2) - neigh_length_c, 1)
                     : min(xyz_vals(rp_ind[jp], 2) + neigh_length_c, N),
                       ... max(xyz_vals(rp_ind[jp], 3) - neigh_length_c, 1)
                     : min(xyz_vals(rp_ind[jp], 3) + neigh_length_c, N));
  cind = ind_box(find(ind_box ~ = 0));
  cind = cind(
      find(full_box_mean(cind) >
           pp_lambda -
               1.96 * sqrt(pp_lambda / std::pow((2 * neigh_length_m + 1), 3))));
  connection_map(cind( :)) = 1;
  end connection_map(rp_ind) = 0;
  connection_map(gp_ind) = 0;
  rp_ind = union(rp_ind, find(connection_map == 1));
  // This last set of Possion Process index includes all points not in total
  // signal //
  gp_rp_ind = union(gp_ind, rp_ind);
  pp_ind = setdiff((1:N_ind)',gp_rp_ind);




  // Determine Noise in Signal //
  std::vector<int> ref_hist;
  for (int j1 = 1; j1 <= static_cast<int>(event_max); j1++) {
    int match;
    for (int j2 = 0; j2 < N_ind; j2++) {
      if (static_cast<int>(event_vals[j2]) == j1) match++;
    }
    ref_hist.push_back(match);
  }
  double ref_hist =
round(mean(full_box_density(pp_ind))*accumarray(event_vals(gp_rp_ind),1)/mean(full_box_density(gp_rp_ind)));
  event_counts.clear();
  for (int j1 = 1; j1 <= static_cast<int>(ref_hist.size(); j1++) {
    event_counts.push_back(j1);
  }
  double ns_dist =
(std::pow((pp_lambda),event_counts))*exp(-pp_lambda)/gsl_sf_fact(event_counts);
  double ns_N = (transpose(ns_dist)*ref_hist)/(transpose(ns_dist)*ns_dist);
  double ns_hist = round(ns_dist.*ns_N);
  double ns_count = sum(ns_hist.*event_counts);*/
  double ns_count = 0.0;
  double gpSum = 0.0;
  double rpSum = 0.0;
  double ppSum = 0.0;
  for (int j1 = 0; j1 < N_ind; j1++) {
    if (find(gp_ind.begin(), gp_ind.end(), j1) != gp_ind.end())
      gpSum += event_vals[j1];
    if (find(rp_ind.begin(), rp_ind.end(), j1) != rp_ind.end())
      rpSum += event_vals[j1];
    if (find(pp_ind.begin(), pp_ind.end(), j1) != pp_ind.end())
      ppSum += event_vals[j1];
  }
  double sum3D = 0.0;
  for (int j1 = 0; j1 < gridPts[0]; j1++) {
    for (int j2 = 0; j2 < gridPts[1]; j2++) {
      for (int j3 = 0; j3 < gridPts[2]; j3++) {
        std::vector<double> tmp;
        tmp.push_back(j1 - c_gp_parm[0]);
        tmp.push_back(j2 - c_gp_parm[1]);
        tmp.push_back(j3 - c_gp_parm[2]);
        tmp = XYZ_est * tmp;
        double gp = exp(-std::pow(tmp[0], 2) / c_gp_parm[6]) *
                    exp(-std::pow(tmp[1], 2) / c_gp_parm[7]) *
                    exp(-std::pow(tmp[2], 2) / c_gp_parm[8]);
        sum3D += best_norm * gp / gp_dist_sum;
      }
    }
  }

  std::cout << best_err << "  " << sum3D << "  " << gpSum << "  " << rpSum
            << "  " << ppSum << "\n";
  double Full_Signal = gpSum;       // + rpSum;
  double Background_in_Signal = 0.; // ns_count;
  double Background_Signal = ppSum;

  intensity = Full_Signal - Background_in_Signal;
  errorSquared = Full_Signal + Background_in_Signal;
  return;
}

/**
 * Runs the BinMD algorithm on the input to provide the output workspace
 * All slicing algorithm properties are passed along
 * @return MDHistoWorkspace as a result of the binning
 */
MDHistoWorkspace_sptr
IntegratePeaks3DFit::binEvent(double Qx, double Qy, double Qz, double box,
                              int gridPts, const IMDWorkspace_sptr &ws) {
  IAlgorithm_sptr binMD = createChildAlgorithm("BinMD", 0.0, 0.3);
  binMD->setProperty("InputWorkspace", ws);
  binMD->setProperty("AlignedDim0",
                     ws->getDimension(0)->getName() + "," +
                         boost::lexical_cast<std::string>(Qx - box) + "," +
                         boost::lexical_cast<std::string>(Qx + box) + "," +
                         std::to_string(gridPts));
  binMD->setProperty("AlignedDim1",
                     ws->getDimension(1)->getName() + "," +
                         boost::lexical_cast<std::string>(Qy - box) + "," +
                         boost::lexical_cast<std::string>(Qy + box) + "," +
                         std::to_string(gridPts));
  binMD->setProperty("AlignedDim2",
                     ws->getDimension(2)->getName() + "," +
                         boost::lexical_cast<std::string>(Qz - box) + "," +
                         boost::lexical_cast<std::string>(Qz + box) + "," +
                         std::to_string(gridPts));
  binMD->setPropertyValue("AxisAligned", "1");
  binMD->setPropertyValue("OutputWorkspace", "out");
  binMD->executeAsChildAlg();
  Workspace_sptr outputWS = binMD->getProperty("OutputWorkspace");
  return boost::dynamic_pointer_cast<MDHistoWorkspace>(outputWS);
}

/**
 * Runs the BinMD algorithm on the input to provide the output workspace
 * All slicing algorithm properties are passed along
 * @return MDHistoWorkspace as a result of the binning
 */
MDHistoWorkspace_sptr
IntegratePeaks3DFit::cropHisto(double Qx, double Qy, double Qz, double box,
                               const IMDWorkspace_sptr &ws) {
  IAlgorithm_sptr cropMD =
      createChildAlgorithm("IntegrateMDHistoWorkspace", 0.0, 0.3);
  cropMD->setProperty("InputWorkspace", ws);

  cropMD->setProperty("P1Bin", boost::lexical_cast<std::string>(Qx - box) +
                                   ",0," +
                                   boost::lexical_cast<std::string>(Qx + box));
  cropMD->setProperty("P2Bin", boost::lexical_cast<std::string>(Qy - box) +
                                   ",0," +
                                   boost::lexical_cast<std::string>(Qy + box));
  cropMD->setProperty("P3Bin", boost::lexical_cast<std::string>(Qz - box) +
                                   ",0," +
                                   boost::lexical_cast<std::string>(Qz + box));

  cropMD->setPropertyValue("OutputWorkspace", "out");
  cropMD->executeAsChildAlg();
  IMDHistoWorkspace_sptr outputWS = cropMD->getProperty("OutputWorkspace");
  return boost::dynamic_pointer_cast<MDHistoWorkspace>(outputWS);
}

} // namespace MDAlgorithms
} // namespace Mantid
