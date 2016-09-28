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
#include "MantidCurveFitting/FortranDefs.h"
#include <boost/math/special_functions/round.hpp>
#include <gsl/gsl_sf_gamma.h> // for factorial
#include <gsl/gsl_linalg.h> // for SVD
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
using namespace CurveFitting;

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(IntegratePeaks3DFit)

//----------------------------------------------------------------------------------------------
/**
  * Initialize the algorithm's properties.
  */
void IntegratePeaks3DFit::init() {
  declareProperty(
      make_unique<WorkspaceProperty<IMDWorkspace>>("InputWorkspace", "",
                                                   Direction::Input),
      "An input Sample MDHistoWorkspace or MDEventWorkspace.");
  declareProperty("Delta", 0.5,
                  "Distance from center to edge of box.");
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
      make_unique<WorkspaceProperty<>>("FluxWorkspace", "", Direction::Input,
                                       PropertyMode::Optional, fluxValidator),
      "An optional input workspace containing momentum dependent flux for "
      "normalization.");
  declareProperty(make_unique<WorkspaceProperty<>>(
                      "SolidAngleWorkspace", "", Direction::Input,
                      PropertyMode::Optional, solidAngleValidator),
                  "An optional input workspace containing momentum integrated "
                  "vanadium for normalization "
                  "(a measure of the solid angle).");

  declareProperty(make_unique<WorkspaceProperty<PeaksWorkspace>>(
                      "PeaksWorkspace", "", Direction::Input),
                  "A PeaksWorkspace containing the peaks to integrate.");

  declareProperty(
      make_unique<WorkspaceProperty<PeaksWorkspace>>("OutputWorkspace", "",
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
  PARALLEL_FOR1(peakWS)
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
      histoBox = normalize(pos.X(), pos.Y(), pos.Z(), box, gridPts, flux, sa, m_eventWS);
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
IntegratePeaks3DFit::normalize(double h, double k, double l, double box, int gridPts,
                               const MatrixWorkspace_sptr &flux,
                               const MatrixWorkspace_sptr &sa,
                               const IMDEventWorkspace_sptr &ws) {
  IAlgorithm_sptr normAlg = createChildAlgorithm("MDNormSCD");
  normAlg->setProperty("InputWorkspace", ws);
  normAlg->setProperty("AlignedDim0", ws->getDimension(0)->getName() +
                       "," + boost::lexical_cast<std::string>(h - box) +
                           "," + boost::lexical_cast<std::string>(h + box) +
                           "," + std::to_string(gridPts));
  normAlg->setProperty("AlignedDim1", ws->getDimension(1)->getName() +
                       "," + boost::lexical_cast<std::string>(k - box) +
                           "," + boost::lexical_cast<std::string>(k + box) +
                           "," + std::to_string(gridPts));
  normAlg->setProperty("AlignedDim2", ws->getDimension(2)->getName() +
                       "," + boost::lexical_cast<std::string>(l - box) +
                           "," + boost::lexical_cast<std::string>(l + box) +
                           "," + std::to_string(gridPts));
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

double eps = 2.2204e-16;           // Floating-point relative accuracy
int do_material = 2;            // Flag to set material files to analyis
                            // do_material = 1 > Silicon
                            // do_material = 2 > Sc
int n_break_its = 100;          // Max number of iterations for each optimization
int max_sigma = 30;             // Max STD in Pixels
int min_sigma = 1;              // Min STD in Pixels
int max_mu = 10;                // Max deviation from Center in Pixels
int N_event_thres = 20;         // Minimum number of events to determine signal


// Begin Program
int neigh_length_c = 0;
int neigh_length_m = 0;

if (do_material == 1) {
    // Silicon file info//
    neigh_length_c = 3;         // Neighborhood size for connection map is 2*neigh_length_c+1
    neigh_length_m = 3;         // Neighborhood size for mean calculation is 2*neigh_length_m+1
}
else if (do_material == 2) {
    // Uncomment for Sc //
    neigh_length_c = 1;         // Neighborhood size for connection map is 2*neigh_length_c+1
    neigh_length_m = 3;         // Neighborhood size for mean calculation is 2*neigh_length_m+1
}
else {
    std::cout <<"Material Flag not defined\n";
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
std::vector<int> x_vals,y_vals,z_vals;
std::vector<double> event_vals, sigma_events;
double event_tot= 0;
double event_max = 0;
int N = gridPts[0];
double N2 = static_cast<double>(N)/2.0;
int N_ind = 0;
for (int j1 = 0; j1 <  gridPts[0]; j1++) {
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
      N_ind++;
      }
    }
  }
}
std::vector<int>event_hist;
std::vector<int>event_counts;
int N_non_zeros = 0;
for (int j1 = 1; j1 <=  static_cast<int>(event_max); j1++) {
  int match;
  for (int j2 = 0; j2 < N_ind; j2++) {
    if(static_cast<int>(event_vals[j2]) == j1) match++;
  }
  event_counts.push_back(j1);
  event_hist.push_back(match);
  if (match > 0) N_non_zeros++;
}
   double pp_lambda = 0.;

    if (N_ind < N_event_thres) {
        std::cout <<"Not enough signal to determine distribution\n";
    }
    else {
        // Initialize Estimate of Poisson Distribution //
        if (N_non_zeros>1) {
            pp_lambda = std::max(eps,(std::pow((event_hist[2]/gsl_sf_fact(event_counts[2]))/(event_hist[1]/gsl_sf_fact(event_counts[1])),(1.0/(event_counts[2]-event_counts[1])))));
        std::vector<double> pp_dist;
        for (int j1 = 0; j1 <  N_ind; j1++) {
            pp_dist.push_back( (std::pow((pp_lambda),event_counts[j1]))*exp(-pp_lambda)/gsl_sf_fact(event_counts[j1]));
        }
            double pp_N = std::inner_product(pp_dist.begin(), pp_dist.end(), event_hist.begin(), 0.0)/std::inner_product(pp_dist.begin(), pp_dist.end(), pp_dist.begin(), 0.0);
            double best_pp_err = 0.;
            for (int j1 = 0; j1 <  N_ind; j1++) {
            best_pp_err += std::pow((pp_N*pp_dist[j1]-event_hist[j1]),2);
            }
            best_pp_err /= event_tot;

            double d_pp_lambda = .1*pp_lambda;
            double did_mv = 0;
            
            for (int jpit = 1; jpit <  n_break_its; jpit++) {
                for (int js  = 1; js <  2; js++) {
                    double c_pp_lambda = pp_lambda;
                    c_pp_lambda += (std::pow((-1),js))*d_pp_lambda;
                    std::vector<double> pp_dist;
                    for (int j1 = 0; j1 <  N_ind; j1++) {
                        pp_dist.push_back( (std::pow((c_pp_lambda),event_counts[j1]))*exp(-c_pp_lambda)/gsl_sf_fact(event_counts[j1]));
                    }
                        double pp_N = std::inner_product(pp_dist.begin(), pp_dist.end(), event_hist.begin(), 0.0)/std::inner_product(pp_dist.begin(), pp_dist.end(), pp_dist.begin(), 0.0);
                    double c_best_pp_err = 0.;
                    for (int j1 = 0; j1 <  N_ind; j1++) {
                    c_best_pp_err += std::pow((pp_N*pp_dist[j1]-event_hist[j1]),2);
                    }
                    c_best_pp_err /= event_tot;
                    if (c_best_pp_err<best_pp_err) {
                        best_pp_err = c_best_pp_err;
                        pp_lambda = c_pp_lambda;
                        did_mv = 1;
                    }
                }
                if (did_mv ==0)
                    d_pp_lambda = d_pp_lambda/2;
                else
                    d_pp_lambda = 1.2*d_pp_lambda;
                if (d_pp_lambda/pp_lambda<1e-4) 
                    break;
            }
        }
      }
        

  std::vector<double> full_box_mean;
  std::vector<double> full_box_density;
  for (int j0 = 0; j0 <  N_ind; j0++) {
    double meanB = 0;
    int num0 = 0;
    int num = 0;
    for (int j1 = x_vals[j0]-neigh_length_m; j1 <= x_vals[j0]+neigh_length_m; j1++) {
      if (j1 < 0 || j1 >= gridPts[0]) continue;
      for (int j2 = y_vals[j0]-neigh_length_m; j2 <= y_vals[j0]+neigh_length_m; j2++) {
        if (j2 < 0 || j2 >= gridPts[1]) continue;
        for (int j3 = z_vals[j0]-neigh_length_m; j3 <= z_vals[j0]+neigh_length_m; j3++) {
          if (j3 < 0 || j3 >= gridPts[2]) continue;
          int iPts = j1 + gridPts[0] * (j2 + gridPts[1] * j3);
          meanB += num_events[iPts];
          num0 ++;
          if(num_events[iPts] > 0.) num ++;
        }
      }
    }
    full_box_mean.push_back(meanB/static_cast<double>(num0));
    full_box_density.push_back(static_cast<double>(num/static_cast<double>(num0)));
  } 
  std::vector<double> conditional_event_vals;
  double conditional_event_vals_tot;
  for (int j0 = 0; j0 <  N_ind; j0++) {
       if(full_box_mean[j0]>pp_lambda+1.96*sqrt(pp_lambda/std::pow((2*neigh_length_m+1),3)) &&
           std::max(fabs(z_vals[j0]-N2),std::max(fabs(x_vals[j0]-N2),fabs(y_vals[j0]-N2)))>max_mu+2*max_sigma) {  // Remove points outside of prescribed domain
         conditional_event_vals.push_back(event_vals[j0]);
         conditional_event_vals_tot += event_vals[j0];
       }
       else conditional_event_vals.push_back(0.0);
  }
        // Estimate initial parameters for Gaussian distribution //
        std::vector<double> mu_est_x;
        std::vector<double> mu_est_y;
        std::vector<double> mu_est_z;
        double muMax= max_mu;
        for (int j0 = 0; j0 <  N_ind; j0++) {
          mu_est_x.push_back(x_vals[j0]*conditional_event_vals[j0]/conditional_event_vals_tot);
          mu_est_y.push_back(y_vals[j0]*conditional_event_vals[j0]/conditional_event_vals_tot);
          mu_est_z.push_back(z_vals[j0]*conditional_event_vals[j0]/conditional_event_vals_tot);
          if (fabs(mu_est_x[j0]-N2 > muMax)) muMax = fabs(mu_est_x[j0]-N2);
          if (fabs(mu_est_y[j0]-N2 > muMax)) muMax = fabs(mu_est_y[j0]-N2);
          if (fabs(mu_est_z[j0]-N2 > muMax)) muMax = fabs(mu_est_z[j0]-N2);
        }
        if (muMax >max_mu) {
          for (int j0 = 0; j0 <  N_ind; j0++) {
            mu_est_x[j0] = 0.5;
            mu_est_y[j0] = 0.5;
            mu_est_z[j0] = 0.5;
          }
        }
        
        DoubleFortranMatrix  covar_mat_est;
        for (int j0 = 0; j0 <  N_ind; j0++) {
            covar_mat_est[0][0] += (conditional_event_vals(jp))*(x_vals(jp)-mu_est_x(jp))*(x_vals(jp)-mu_est_x(jp));
            covar_mat_est[1][0] += (conditional_event_vals(jp))*(x_vals(jp)-mu_est_y(jp))*(y_vals(jp)-mu_est_x(jp));
            covar_mat_est[2][0] += (conditional_event_vals(jp))*(x_vals(jp)-mu_est_z(jp))*(z_vals(jp)-mu_est_x(jp));
            covar_mat_est[0][1] += (conditional_event_vals(jp))*(y_vals(jp)-mu_est_x(jp))*(x_vals(jp)-mu_est_y(jp));
            covar_mat_est[1][1] += (conditional_event_vals(jp))*(y_vals(jp)-mu_est_y(jp))*(y_vals(jp)-mu_est_y(jp));
            covar_mat_est[2][1] += (conditional_event_vals(jp))*(y_vals(jp)-mu_est_z(jp))*(z_vals(jp)-mu_est_y(jp));
            covar_mat_est[0][2] += (conditional_event_vals(jp))*(z_vals(jp)-mu_est_x(jp))*(x_vals(jp)-mu_est_z(jp));
            covar_mat_est[1][2] += (conditional_event_vals(jp))*(z_vals(jp)-mu_est_y(jp))*(y_vals(jp)-mu_est_z(jp));
            covar_mat_est[2][2] += (conditional_event_vals(jp))*(z_vals(jp)-mu_est_z(jp))*(z_vals(jp)-mu_est_z(jp));
        }
        for (int j0 = 0; j0 <  3; j0++) {
          for (int j1 = 0; j1 <  3; j1++)
          covar_mat_est[j0][j1] = covar_mat_est[j0][j1]/conditional_event_vals_tot;
        }
        auto n = 3;
        DoubleFortranMatrix u_est = covar_mat_est;
        DoubleFortranMatrix v_est(n, n);
        DoubleFortranVector d_est(n);
        DoubleFortranVector work(n);
        gsl_linalg_SV_decomp(u_est.gsl(), v_est.gsl(), d_est.gsl(), work.gsl());
        
        theta_x_est = atan2(u_est(3,2),u_est(3,3));
        theta_y_est = atan2(-u_est(3,1),sqrt(std::pow(u_est(3,2),2)+std::pow(u_est(3,3),2)));
        theta_z_est = atan2(u_est(2,1),u_est(1,1));
        sigma = d_est;
        
       /* for (int jp = 0; jp <  3; jp++) {
            if sqrt(sigma(jp)/2)<min_sigma
                sigma(jp) = 2.1*std::pow((min_sigma),2);
            elseif sqrt(sigma(jp)/2)>max_sigma
                sigma(jp) = 1.9*std::pow((max_sigma),2);
            end
        end

        gp_parm = [mu_est(:);[theta_x_est;theta_y_est;theta_z_est];sigma];
        c_gp_parm = gp_parm;
        
        X_est = [[1 0 0]; [0 cos(c_gp_parm(4)) -sin(c_gp_parm(4))]; [0 sin(c_gp_parm(4)) cos(c_gp_parm(4))]];
        Y_est = [[cos(c_gp_parm(5)) 0 sin(c_gp_parm(5))]; [0 1 0]; [-sin(c_gp_parm(5)) 0 cos(c_gp_parm(5))]];
        Z_est = [[cos(c_gp_parm(6)) -sin(c_gp_parm(6)) 0]; [sin(c_gp_parm(6)) cos(c_gp_parm(6)) 0]; [0 0 1]];
        rot_cord = (xyz_vals-repmat(c_gp_parm(1:3)',N_ind,1))*Z_est*Y_est*X_est;
        
        
        f_x = exp(std::pow(-rot_cord(:,1),2)/c_gp_parm(7));
        f_y = exp(std::pow(-rot_cord(:,2),2)/c_gp_parm(8));
        f_z = exp(std::pow(-rot_cord(:,3),2)/c_gp_parm(9));
        
        gp_dist = f_x.*f_y.*f_z;
        gp_dist = gp_dist/sum(gp_dist(:));
        best_norm = (transpose(gp_dist)*conditional_event_vals)/(transpose(gp_dist)*gp_dist);
        gp_events = round(best_norm*gp_dist);
        best_err = sqrt(sum(std::pow((conditional_event_vals-best_norm*gp_dist),2)))/event_tot;
        
        // Estimate Gaussian distribution given Poisson distribution and Extreme Value Distribution //
        hfg = figure;
        trackbest_err = best_err;
        d_gp_parm = .01*fabs(gp_parm);
        d_gp_parm = max(d_gp_parm,.1);      // Minimal starting search distance of optimization parameters.
        
        for jpit = 1:n_break_its
            did_mv_vec = zeros(length(gp_parm),1);
            for jp = 1:length(gp_parm)
                for js = 1:2
                    c_gp_parm = gp_parm;
                    c_gp_parm(jp) = c_gp_parm(jp)+(std::pow((-1),js))*d_gp_parm(jp);
                    do_opt = 1;
                    // Enforce paramameter maximum domains %
                    if jp<4
                        if fabs(c_gp_parm(jp)-N2)>max_mu
                            do_opt = 0;
                        end
                    end
                    if jp>6
                        if (sqrt(c_gp_parm(jp)/2)<min_sigma | sqrt(c_gp_parm(jp)/2)>max_sigma)
                            do_opt = 0;
                        end
                    end
                    
                    if do_opt == 1
                        X_est = [[1 0 0]; [0 cos(c_gp_parm(4)) -sin(c_gp_parm(4))]; [0 sin(c_gp_parm(4)) cos(c_gp_parm(4))]];
                        Y_est = [[cos(c_gp_parm(5)) 0 sin(c_gp_parm(5))]; [0 1 0]; [-sin(c_gp_parm(5)) 0 cos(c_gp_parm(5))]];
                        Z_est = [[cos(c_gp_parm(6)) -sin(c_gp_parm(6)) 0]; [sin(c_gp_parm(6)) cos(c_gp_parm(6)) 0]; [0 0 1]];
                        rot_cord = (xyz_vals-repmat(c_gp_parm(1:3)',N_ind,1))*Z_est*Y_est*X_est;
                        
                        f_x = exp(std::pow(-rot_cord(:,1),2)/c_gp_parm(7));
                        f_y = exp(std::pow(-rot_cord(:,2),2)/c_gp_parm(8));
                        f_z = exp(std::pow(-rot_cord(:,3),2)/c_gp_parm(9));
                        c_gp_dist = f_x.*f_y.*f_z;
                        c_gp_dist = c_gp_dist/sum(c_gp_dist(:));
                        best_norm = (transpose(c_gp_dist)*conditional_event_vals)/(transpose(c_gp_dist)*c_gp_dist);
                        c_gp_events = round(best_norm.*c_gp_dist);
                        c_best_err = sqrt(sum(std::pow((conditional_event_vals-best_norm*c_gp_dist),2)))/event_tot;
                        if c_best_err<best_err
                            best_err = c_best_err;
                            gp_parm = c_gp_parm;
                            gp_dist = c_gp_dist;
                            gp_events = c_gp_events;
                            did_mv_vec(jp) = 1;
                        end
                    end
                end
                if did_mv_vec(jp) ==0
                    d_gp_parm(jp) = d_gp_parm(jp)/2;
                else
                    d_gp_parm(jp) = 1.2*d_gp_parm(jp);
                end
            end
            
            figure(hfg);
            trackbest_err = [trackbest_err;best_err];
            semilogy(trackbest_err)
            set(gca,'FontSize',16);
            //title(['Error for Multivariate (',num2str(best_err,'%10.3e\n'),') -- Gaussian, Slice ',num2str(jtarg)],'FontSize',16)
            drawnow
            
            if max(d_gp_parm/fabs(gp_parm))<1e-3
                break
            end
        end
        
        close(hfg);
        
        // Estimate Poisson distribution given Gaussian distribution //
        
        conditional_event_vals = max(event_vals-gp_events,0);
        pp_ind = find(round(gp_events)<1);
        gp_ind = find(round(gp_events)>0);
        event_hist = accumarray(event_vals(pp_ind),1);
        event_counts = (1:length(event_hist))';
        
        N_non_zeros = length(find(event_hist~=0));
        if N_non_zeros>1
            pp_lambda = max(eps,std::pow(((event_hist[2]/gsl_sf_fact(event_counts[2]))/(event_hist[1]/gsl_sf_fact(event_counts[1]))),(1/(event_counts[2]-event_counts[1]))));
            pp_dist = (std::pow((pp_lambda),event_counts))*exp(-pp_lambda)/gsl_sf_fact(event_counts);
            pp_N = (transpose(pp_dist)*event_hist)/(transpose(pp_dist)*pp_dist);
            best_pp_err = sum(std::pow((pp_N*pp_dist-event_hist),2))/event_tot;
            d_pp_lambda = .1*pp_lambda;
            
            for jpit = 1:n_break_its
                did_mv = 0;
                for js = 1:2
                    c_pp_lambda = pp_lambda;
                    c_pp_lambda = c_pp_lambda+(std::pow((-1),js))*d_pp_lambda;
                    pp_dist = (std::pow((c_pp_lambda),event_counts))*exp(-c_pp_lambda)/gsl_sf_fact(event_counts);
                    pp_N = (transpose(pp_dist)*event_hist)/(transpose(pp_dist)*pp_dist);
                    c_best_pp_err = sum(std::pow((pp_N*pp_dist-event_hist,2))/event_tot;
                    if c_best_pp_err<best_pp_err
                        best_pp_err = c_best_pp_err;
                        pp_lambda = c_pp_lambda;
                        did_mv = 1;
                    end
                end
                if did_mv ==0
                    d_pp_lambda = d_pp_lambda/2;
                else
                    d_pp_lambda = 1.2*d_pp_lambda;
                end
                if d_pp_lambda/pp_lambda<1e-4
                    break
                end
            end
            pp_dist = (std::pow((pp_lambda),event_counts))*exp(-pp_lambda)/gsl_sf_fact(event_counts);
            pp_N = (transpose(pp_dist)*event_hist)/(transpose(pp_dist)*pp_dist);
            
            rp_ind = intersect(find(full_box_mean>pp_lambda+1.96*sqrt(pp_lambda/std::pow((2*neigh_length_m+1),3))),pp_ind);
        else
            // Not Enough Information to determine Possion Distribution //
            // Estimate Sparsity of Distribution and use for background //
            pp_lambda = 0;
            rp_ind = [];
        end
        

        // Estimate Residual Distribution given Poisson distribution and Gaussian distribution //

        connection_map = 1:N_ind;
        connection_map(gp_ind) = 0;

        for jp = 1:length(rp_ind);
            ind_box = num_events_ref(max(xyz_vals(rp_ind(jp),1)-neigh_length_c,1):min(xyz_vals(rp_ind(jp),1)+neigh_length_c,N),...
                max(xyz_vals(rp_ind(jp),2)-neigh_length_c,1):min(xyz_vals(rp_ind(jp),2)+neigh_length_c,N),...
                max(xyz_vals(rp_ind(jp),3)-neigh_length_c,1):min(xyz_vals(rp_ind(jp),3)+neigh_length_c,N));
            cind = ind_box(find(ind_box~=0));
            cind = cind(find(full_box_mean(cind)>pp_lambda+1.96*sqrt(pp_lambda/std::pow((2*neigh_length_m+1),3))));
            
            connection_map(cind(:)) = min(connection_map(cind(:)));
        end
        
        // Set of points only in the Tail of distribution.  //
        rp_ind = setdiff(find(connection_map==0),gp_ind);
        
        // Add points at 5// level of pd //
        neig_ind = [];
        connection_map = zeros(N_ind,1);
        for jp = 1:length(rp_ind);
            ind_box = num_events_ref(max(xyz_vals(rp_ind(jp),1)-neigh_length_c,1):min(xyz_vals(rp_ind(jp),1)+neigh_length_c,N),...
                max(xyz_vals(rp_ind(jp),2)-neigh_length_c,1):min(xyz_vals(rp_ind(jp),2)+neigh_length_c,N),...
                max(xyz_vals(rp_ind(jp),3)-neigh_length_c,1):min(xyz_vals(rp_ind(jp),3)+neigh_length_c,N));
            cind = ind_box(find(ind_box~=0));
            cind = cind(find(full_box_mean(cind)>pp_lambda-1.96*sqrt(pp_lambda/std::pow((2*neigh_length_m+1),3))));
            connection_map(cind(:)) = 1;
        end
        connection_map(rp_ind) = 0;
        connection_map(gp_ind) = 0;
        rp_ind = union(rp_ind,find(connection_map==1));
        // This last set of Possion Process index includes all points not in total signal //
        gp_rp_ind = union(gp_ind,rp_ind);
        pp_ind = setdiff((1:N_ind)',gp_rp_ind);
        
        
        
        
        // Determine Noise in Signal //
        ref_hist = round(mean(full_box_density(pp_ind))*accumarray(event_vals(gp_rp_ind),1)/mean(full_box_density(gp_rp_ind)));
        event_counts = (1:length(ref_hist))';
        ns_dist = (std::pow((pp_lambda),event_counts))*exp(-pp_lambda)/gsl_sf_fact(event_counts);
        ns_N = (transpose(ns_dist)*ref_hist)/(transpose(ns_dist)*ns_dist);
        ns_hist = round(ns_dist.*ns_N);
        ns_count = sum(ns_hist.*event_counts);
        
        Full_Signal = sum(event_vals(gp_ind))+sum(event_vals(rp_ind));
        Background_in_Signal = ns_count;
        Background_Signal = sum(event_vals(pp_ind));
        
        //////////////////////////////
        //// Code Below is only for //%
        ////       vis purpose      //%
        //////////////////////////////
        gp_events = zeros(N_ind,1);
        gp_events(gp_ind) = event_vals(gp_ind);
        pp_events = zeros(N_ind,1);
        pp_events(pp_ind) = event_vals(pp_ind);
        rp_events = zeros(N_ind,1);
        rp_events(rp_ind) = event_vals(rp_ind);
        
        ns_events = zeros(N_ind,1);
        ns_inds = [];
        for jp = length(ns_hist):-1:1
            cind = setdiff(find(gp_events+rp_events-ns_events >= jp),ns_inds);
            cind = cind(randperm(length(cind)));
            cind = cind(1:min(ns_hist(jp),length(cind)));
            ns_inds = union(ns_inds,cind);
            ns_events(cind) = jp;
        end
        
        
        
        close all
        h_fig = figure('Position', [0, 0, 1400, 900]);
        cmax = max(num_events(:));
        vis_event(event_vals,num_events,use_ind,Q1,Q2,Q3,h_fig,1,cmax)
        title(['Box, Total Count: ',num2str(sum(event_vals))],'FontSize',12);
        
        vis_event(gp_events(:),num_events,use_ind,Q1,Q2,Q3,h_fig,2,cmax)
        title(['Gaussian Process Signal, Count: ',num2str(sum(gp_events(:)))],'FontSize',12);
        
        vis_event(gp_events+rp_events-ns_events,num_events,use_ind,Q1,Q2,Q3,h_fig,3,cmax)
        title(['Full Signal, Count: ',num2str(sum(gp_events(:)+rp_events(:)-ns_events(:)))],'FontSize',12);
        
        vis_event(pp_events+ns_events,num_events,use_ind,Q1,Q2,Q3,h_fig,4,cmax)
        if do_ref == 1
            title(['Background, Count: ',num2str(sum(pp_events(:)+ns_events(:))),'; BG in Sig: ',num2str(sum(ns_events(:))), '; Reference Signal ',num2str(intensities(jtarg+1)),' ; New Alg ',num2str(sum(gp_events(:)+rp_events(:)-ns_events(:)))],'FontSize',12);
        else
            title(['Background, Count: ',num2str(sum(pp_events(:)+ns_events(:))),'; BG in Sig: ',num2str(sum(ns_events(:)))],'FontSize',12);
        end

        drawnow
//         p=mtit(['Box ',num2str(jtarg),', Signal Decomposition'],'FontSize',14);
//         set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 30])
//         print(h_fig,'-depsc', ['EVE',num2str(jtarg),'.eps'], '-r300');
    end
end



  double minIntensity = Fmin + 0.01 * (Fmax - Fmin);
  int measuredPoints = 0;
  int peakPoints = 0;
  double peakSum = 0.0;
  double measuredSum = 0.0;
  double errSqSum = 0.0;
  double measuredErrSqSum = 0.0;
  for (int Hindex = 0; Hindex < gridPts[0]; Hindex++) {
    for (int Kindex = 0; Kindex < gridPts[1]; Kindex++) {
      for (int Lindex = 0; Lindex < gridPts[2]; Lindex++) {
        int iPts = Hindex + gridPts[0] * (Kindex + gridPts[1] * Lindex);
        if (std::isfinite(F[iPts])) {
          measuredPoints = measuredPoints + 1;
          measuredSum = measuredSum + F[iPts];
          measuredErrSqSum = measuredErrSqSum + SqError[iPts];
          if (F[iPts] > minIntensity) {
            int neighborPoints = 0;
            for (int Hj = -2; Hj < 3; Hj++) {
              for (int Kj = -2; Kj < 3; Kj++) {
                for (int Lj = -2; Lj < 3; Lj++) {
                  int jPts =
                      Hindex + Hj +
                      gridPts[0] * (Kindex + Kj + gridPts[1] * (Lindex + Lj));
                  if (Lindex + Lj >= 0 && Lindex + Lj < gridPts[2] &&
                      Kindex + Kj >= 0 && Kindex + Kj < gridPts[1] &&
                      Hindex + Hj >= 0 && Hindex + Hj < gridPts[0] &&
                      F[jPts] > minIntensity) {
                    neighborPoints = neighborPoints + 1;
                  }
                }
              }
            }
            if (neighborPoints >= neighborPts) {
              peakPoints = peakPoints + 1;
              peakSum = peakSum + F[iPts];
              errSqSum = errSqSum + SqError[iPts];
            }
          }
        } else {
          double minR =
              sqrt(std::pow(float(Hindex) / float(gridPts[0]) - 0.5, 2) +
                   std::pow(float(Kindex) / float(gridPts[1]) - 0.5, 2) +
                   std::pow(float(Lindex) / float(gridPts[0]) - 0.5, 2));
          if (minR < 0.05) {
            intensity = 0.0;
            errorSquared = 0.0;
            return;
          }
        }
      }
    }
  }
  double ratio = float(peakPoints) / float(measuredPoints - peakPoints);
  intensity = peakSum - ratio * (measuredSum - peakSum);
  errorSquared = errSqSum + ratio * ratio * (measuredErrSqSum - errSqSum);
  return;*/
}

/**
 * Runs the BinMD algorithm on the input to provide the output workspace
 * All slicing algorithm properties are passed along
 * @return MDHistoWorkspace as a result of the binning
 */
MDHistoWorkspace_sptr
IntegratePeaks3DFit::binEvent(double Qx, double Qy, double Qz, double box, int gridPts,
                              const IMDWorkspace_sptr &ws) {
  IAlgorithm_sptr binMD = createChildAlgorithm("BinMD", 0.0, 0.3);
  binMD->setProperty("InputWorkspace", ws);
  binMD->setProperty("AlignedDim0", ws->getDimension(0)->getName() +
                     "," + boost::lexical_cast<std::string>(Qx - box) +
                         "," + boost::lexical_cast<std::string>(Qx + box) + "," +
                         std::to_string(gridPts));
  binMD->setProperty("AlignedDim1", ws->getDimension(1)->getName() +
                     "," + boost::lexical_cast<std::string>(Qy - box) +
                         "," + boost::lexical_cast<std::string>(Qy + box) + "," +
                         std::to_string(gridPts));
  binMD->setProperty("AlignedDim2", ws->getDimension(2)->getName() +
                     "," + boost::lexical_cast<std::string>(Qz - box) +
                         "," + boost::lexical_cast<std::string>(Qz + box) + "," +
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

