//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidAlgorithms/Qhelper.h"
#include "MantidKernel/ArrayProperty.h"
#include "MantidKernel/RebinParamsValidator.h"
#include "MantidKernel/UnitFactory.h"
#include "MantidKernel/PhysicalConstants.h"
#include "MantidKernel/VectorHelper.h"
#include "MantidAPI/WorkspaceValidators.h"
#include "MantidAPI/SpectraDetectorMap.h"
#include "MantidDataObjects/Histogram1D.h"

namespace Mantid
{
namespace Algorithms
{


using namespace Kernel;
using namespace API;
using namespace Geometry;


/** Checks if workspaces input to Q1D or Qxy are reasonable
  @param dataWS data workspace
  @param binWS (WavelengthAdj) workpace that will be checked to see if it has one spectrum and the same number of bins as dataWS
  @param detectWS (PixelAdj) passing NULL for this wont raise an error, if set it will be checked this workspace has as many histograms as dataWS each with one bin
  @throw invalid_argument if the workspaces are not mututially compatible
*/
void Qhelper::examineInput(API::MatrixWorkspace_const_sptr dataWS, 
     API::MatrixWorkspace_const_sptr binAdj, API::MatrixWorkspace_const_sptr detectAdj)
{
  if ( dataWS->getNumberHistograms() < 1 )
  {
    throw std::invalid_argument("Empty data workspace passed, can not continue");
  }

  //it is not an error for these workspaces not to exist
  if (binAdj)
  {
    if ( binAdj->getNumberHistograms() != 1 )
    {
      throw std::invalid_argument("The WavelengthAdj workspace must have one spectrum");
    }
    if ( binAdj->readY(0).size() != dataWS->readY(0).size() )
    {
      throw std::invalid_argument("The WavelengthAdj workspace's bins must match those of the detector bank workspace");
    }
    MantidVec::const_iterator reqX = dataWS->readX(0).begin();
    MantidVec::const_iterator testX = binAdj->readX(0).begin();
    for ( ; reqX != dataWS->readX(0).end(); ++reqX, ++testX)
    {
      if ( *reqX != *testX )
      {
        throw std::invalid_argument("The WavelengthAdj workspace must have matching bins with the detector bank workspace");
      }
    }
    if ( binAdj->isDistribution() != dataWS->isDistribution() )
    {
      throw std::invalid_argument("The distrbution/raw counts status of the wavelengthAdj and DetBankWorkspace must be the same, use ConvertToDistribution");
    }
  }
  else if( ! dataWS->isDistribution() )
  {
    //throw std::invalid_argument("The data workspace must be a distrbution if there is no Wavelength dependent adjustment");
  }
  
  if (detectAdj)
  {
    if ( detectAdj->blocksize() != 1 )
    {
      throw std::invalid_argument("The PixelAdj workspace must point to a workspace with single bin spectra, as only the first bin is used");
    }
    if ( detectAdj->getNumberHistograms() != dataWS->getNumberHistograms() )
    {
      throw std::invalid_argument("The PixelAdj workspace must have one spectrum for each spectrum in the detector bank workspace");
    }
  }
}

/** Finds the first index number of the first wavelength bin that should included based on the
*  the calculation: W = Wcut (Rcut-R)/Rcut
*  @param dataWS data workspace
*  @param RCut the radius cut off, should be value of the property RadiusCut
*  @param WCut this wavelength cut off, should be equal to the value WaveCut
*  @param specInd spectrum that is being analysed
*  @return index number of the first bin to include in the calculation
*/
size_t Qhelper::waveLengthCutOff(API::MatrixWorkspace_const_sptr dataWS, const double RCut, const double WCut, 
                                 const size_t specInd) const
{
  double l_WCutOver = 0.0;
  double l_RCut = 0.0;
  if ( RCut > 0 && WCut > 0 )
  {
    l_WCutOver = WCut/RCut;
    l_RCut = RCut;
  }

  if ( !(l_RCut > 0) )
  {
    return 0;
  }
  //get the distance of between this detector and the origin, which should be the along the beam center
  const V3D posOnBank = dataWS->getDetector(specInd)->getPos();
  double R = (posOnBank.X()*posOnBank.X())+(posOnBank.Y()*posOnBank.Y());
  R = std::sqrt(R);

  const double WMin = l_WCutOver*(l_RCut-R);
  const MantidVec & Xs = dataWS->readX(specInd);
  return std::lower_bound(Xs.begin(), Xs.end(), WMin) - Xs.begin();
}

} // namespace Algorithms
} // namespace Mantid

