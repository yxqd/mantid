//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidAlgorithms/GetDetectorOffsets.h"
#include "MantidAPI/WorkspaceValidators.h"
#include "MantidAPI/SpectraDetectorMap.h"
#include "MantidAPI/FileProperty.h"
#include "MantidAPI/FunctionFactory.h"
#include "MantidAPI/IFunction.h"
#include <boost/math/special_functions/fpclassify.hpp>
#include <fstream>
#include <ostream>
#include <iomanip>
#include <sstream>

namespace Mantid
{
  namespace Algorithms
  {

    // Register the class into the algorithm factory
    DECLARE_ALGORITHM(GetDetectorOffsets)

    using namespace Kernel;
    using namespace API;

    /// Constructor
    GetDetectorOffsets::GetDetectorOffsets() :
      API::Algorithm()
    {}

    /// Destructor
    GetDetectorOffsets::~GetDetectorOffsets()
    {}

    /** Initialisation method. Declares properties to be used in algorithm.
     *
     */
    void GetDetectorOffsets::init()
    {

      declareProperty(new WorkspaceProperty<>("InputWorkspace","",Direction::Input,
          new WorkspaceUnitValidator<>("dSpacing")),"A 2D workspace with X values of d-spacing");
      declareProperty(new API::WorkspaceProperty<>("OutputWorkspace","",Direction::Output),"Workspace containing the offsets");
      BoundedValidator<double> *mustBePositive = new BoundedValidator<double>();
      mustBePositive->setLower(0);

      declareProperty("Step",0.001, mustBePositive,
        "Step size used to bin d-spacing data");
      declareProperty("DReference",2.0, mustBePositive->clone(),
         "Center of reference peak in d-space");
      declareProperty("XMin",0.0, "Minimum of CrossCorrelation data to search for peak, usually negative");
      declareProperty("XMax",0.0, "Maximum of CrossCorrelation data to search for peak, usually positive");
      declareProperty(new FileProperty("GroupingFileName","", FileProperty::Save, ".cal"),
		      "The name of the output CalFile" );
    }

    /** Executes the algorithm
     *
     *  @throw Exception::FileError If the grouping file cannot be opened or read successfully
     */
    void GetDetectorOffsets::exec()
    {
      retrieveProperties();
      // Fit all the spectra with a gaussian
      std::string filename=getProperty("GroupingFileName");
      const SpectraDetectorMap& specMap = inputW->spectraMap();
      std::fstream out;
      out.open(filename.c_str(),std::ios::out);
      if (!out.is_open())
      {
        std::runtime_error("Problem opening file"+filename);
      }
      else
      {
        int n=0;
        out << "# Offsets generated by Mantid" << std::endl;
        for (int i=0;i<nspec;++i)
        {
          outputW->getAxis(1)->spectraNo(i)=inputW->getAxis(1)->spectraNo(i);
          double offset=fitSpectra(i);
          // Assign the value of the offset
          outputW->dataY(i)[0]=offset;
          // Put it into file
          int specno=outputW->getAxis(1)->spectraNo(i);
          const std::vector<int> dets=specMap.getDetectors(specno);
          for (unsigned int j=0;j<dets.size();j++)
          {
            out << std::fixed << std::setw(9) << n++ <<
                std::fixed << std::setw(15) << dets[j] <<
                std::fixed << std::setprecision(7) << std::setw(15) << offset <<
                std::fixed << std::setw(8) << "1" <<
                std::fixed << std::setw(8) << "1"  << "\n";
          }
          progress((double)(i)/nspec);
        }
        out.close();
      }

      setProperty("OutputWorkspace",outputW);
    }

    void GetDetectorOffsets::retrieveProperties()
    {
      inputW=getProperty("InputWorkspace");
      Xmin=getProperty("XMin");
      Xmax=getProperty("XMax");
      if (Xmin>=Xmax)
        throw std::runtime_error("Must specify Xmin<Xmax");
      dreference=getProperty("DReference");
      step=getProperty("Step");
      nspec=inputW->getNumberHistograms();
      outputW=API::WorkspaceFactory::Instance().create(inputW,nspec,2,1);
    }

   /** Calls Gaussian1D as a child algorithm to fit the offset peak in a spectrum
    *  @param s The spectrum index to fit
    *  @return The calculated offset value
    */
    double GetDetectorOffsets::fitSpectra(const int s)
    {
      // Find point of peak centre
      const MantidVec & yValues = inputW->readY(s);
      MantidVec::const_iterator it = std::max_element(yValues.begin(), yValues.end());
      const double peakHeight = *it; 
      const double peakLoc = inputW->readX(s)[it - yValues.begin()];
      // Return offset of 0 if peak of Cross Correlation is nan (Happens when spectra is zero)
      if ( boost::math::isnan(peakHeight) ) return (0.);

      IAlgorithm_sptr fit_alg;
      try
      {
        //set the subalgorithm no to log as this will be run once per spectra
        fit_alg = createSubAlgorithm("Fit",-1,-1,false);
      } catch (Exception::NotFoundError&)
      {
        g_log.error("Can't locate Gaussian1D");
        throw ;
      }
      fit_alg->setProperty("InputWorkspace",inputW);
      fit_alg->setProperty("WorkspaceIndex",s);
      fit_alg->setProperty("StartX",Xmin);
      fit_alg->setProperty("EndX",Xmax);
      fit_alg->setProperty("MaxIterations",100);

      std::ostringstream fun_str;
      fun_str << "name=LinearBackground;name=Gaussian,Height="<<peakHeight<<",";
      fun_str << "PeakCentre="<<peakLoc<<",Sigma=10.0";

      fit_alg->setProperty("Function",fun_str.str());

      try
      {
        fit_alg->execute();
      }
      catch (std::runtime_error&)
      {
        g_log.error("Unable to successfully run Gaussian1D sub-algorithm");
        throw;
      }

      if ( ! fit_alg->isExecuted() )
      {
        g_log.error("Unable to successfully run Gaussian1D sub-algorithm");
        throw std::runtime_error("Unable to successfully run Gaussian1D sub-algorithm");
      }

      IFunction* fun = FunctionFactory::Instance().createInitialized(fit_alg->getPropertyValue("Function"));
      if (!fun)
      {
        throw std::runtime_error("FunctionFactory cannot create function returned from Fit algorithm");
      }
      const double offset = fun->getParameter("f1.PeakCentre");
      delete fun;
      return (-offset*step/(dreference+offset*step));
    }



  } // namespace Algorithm
} // namespace Mantid
