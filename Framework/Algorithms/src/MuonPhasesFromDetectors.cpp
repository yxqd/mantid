#include "MantidAlgorithms/MuonPhasesFromDetectors.h"

#include "MantidAPI/Axis.h"
#include "MantidAPI/FunctionFactory.h"
#include "MantidAPI/GroupingLoader.h"
#include "MantidAPI/ITableWorkspace.h"
#include "MantidAPI/IFunction.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/MultiDomainFunction.h"
#include "MantidAPI/Run.h"
#include "MantidAPI/TableRow.h"
#include "MantidIndexing/IndexInfo.h"
#include "MantidKernel/ArrayProperty.h"
#include "MantidKernel/PhysicalConstants.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/WorkspaceGroup.h"

#include "MantidKernel/ListValidator.h"
#include <math.h>



namespace Mantid {
namespace Algorithms {

using namespace Kernel;
namespace {
	// create instrument workspace
	API::MatrixWorkspace_sptr createInstrumentWS(API::IAlgorithm_sptr instrumentWS,std::string instrument){
		instrumentWS->initialize();
		instrumentWS->setProperty("Instrument", instrument);
		instrumentWS->setProperty("BinParams", "0,1,32");
		instrumentWS->setProperty("OutputWorkspace", "tmp");
		instrumentWS->execute();
		return instrumentWS->getProperty("OutputWorkspace");
	}
	// creates an empty phase table
	API::ITableWorkspace_sptr createPhaseTable(API::IAlgorithm_sptr tableAlg) {
		tableAlg->initialize();
		tableAlg->setProperty("OutputWorkspace", "table");
		tableAlg->execute();
		//add headings to table
		API::ITableWorkspace_sptr phaseTable = tableAlg->getProperty("OutputWorkspace");
		phaseTable->addColumn("int", "DetectorID");
		phaseTable->addColumn("double", "Phase");
		phaseTable->addColumn("double", "Asymmetry");
		return phaseTable;
	}

	// Now convert phases to interval [0, 2PI)
	double getPhaseFromZeroTo2Pi(const double phi) {
		int factor = static_cast<int>(floor(phi / (2. * M_PI)));
		if (factor) {
			return phi - factor * 2. * M_PI;
		}
		return phi;
	}
	// Return asymmetry and phase pair for field along x axis
	std::pair<double, double> xField(const V3D det, const double radius) {
		double phi = atan2(det.Z(), det.Y()) - M_PI;
		double asymm = sqrt(det.Z()*det.Z() + det.Y()*det.Y()) / radius;
		std::pair <double, double> asymmPhasePair(asymm, getPhaseFromZeroTo2Pi(phi));
		return asymmPhasePair;
	}
	// Return asymmetry and phase pair for field along y axis
	std::pair<double, double> yField(const V3D det, const double radius) {
		double phi = atan2(det.X(), det.Z()) - M_PI;
		double asymm = sqrt(det.X()*det.X() + det.Z()*det.Z()) / radius;
		std::pair <double, double> asymmPhasePair(asymm, getPhaseFromZeroTo2Pi(phi));
		return asymmPhasePair;
	}
	// Return asymmetry and phase pair for field along z axis
	std::pair<double, double> zField(const V3D det, const double radius) {
		double phi = atan2(det.Y(), det.X())-M_PI;
		double asymm = sqrt(det.Y()*det.Y() + det.X()*det.X()) / radius;
		std::pair <double, double> asymmPhasePair(asymm, getPhaseFromZeroTo2Pi(phi));
		return asymmPhasePair;
	}
	// returns a pointer to relevant function to calculat asymmetry, pahse pair
	typedef std::pair<double,double> (*fptr)(const V3D det, const double radius);
	fptr getAsymmPhasePair(std::string axis) {
		if (axis == "X") {
			return xField;
		}
		else if (axis == "Y") {
			return yField;
		}
		// Z
		else {
			return zField;
		}
	}

}
// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(MuonPhasesFromDetectors)

/** Initializes the algorithm's properties.
 */
void MuonPhasesFromDetectors::init() {
	Kernel::IValidator_sptr instrumentList = boost::make_shared<Mantid::Kernel::StringListValidator>(
		std::vector<std::string>{"MUSR", "EMU", "HIFI", "ARGUS", "CHRONUS", "MUT"});
  declareProperty("Instrument", "MUSR",instrumentList,
	  "The instrument to calculate phases for");  

  Kernel::IValidator_sptr field=boost::make_shared<Mantid::Kernel::StringListValidator>(
	  std::vector<std::string>{"X", "Y", "Z"});
  declareProperty("FieldDirection", "X",field,		 
		  "The filed direction");
  declareProperty(make_unique<API::WorkspaceProperty<API::ITableWorkspace>>(
                      "DetectorTable", "", Direction::Output),
                  "Name of the TableWorkspace in which to store the list "
                  "of phases and asymmetries");

}

/** Executes the algorithm.
 */
void MuonPhasesFromDetectors::exec() {
	// create instrument workspace
	API::IAlgorithm_sptr instrumentWS = createChildAlgorithm("CreateSimulationWorkspace");
	std::string instrument = getProperty("Instrument");
	API::MatrixWorkspace_sptr outputWS = createInstrumentWS(instrumentWS, instrument);
	
	// create an empty table
	API::IAlgorithm_sptr tableAlg = createChildAlgorithm("CreateEmptyTableWorkspace");
	API::ITableWorkspace_sptr phaseTable = createPhaseTable(tableAlg);

	std::string axis = getProperty("FieldDirection");
	fptr asymmPhasePairFunction = getAsymmPhasePair(axis);
	int numberHistograms = 0;
	try {
		numberHistograms = boost::numeric_cast<int>(outputWS->getNumberHistograms());
	}catch(boost::numeric::bad_numeric_cast &error) {
		g_log.error("Bad numeric cast for detector ID "+static_cast<std::string>(error.what()));
	}
	for (int j = 0; j < numberHistograms; j++) {
		reportProgress(j, numberHistograms);
		auto det = outputWS->getDetector(j)->getPos() - outputWS->getInstrument()->getSample()->getPos();
		double radius = sqrt(det.X()*det.X() + det.Y()*det.Y() + det.Z()*det.Z());
		std::pair<double, double> aymmPhasePairValues = asymmPhasePairFunction(det, radius);
		API::TableRow row = phaseTable->appendRow();
		row << j << aymmPhasePairValues.second << aymmPhasePairValues.first;
	}

	setProperty("DetectorTable",phaseTable);

}

/**
 * Updates the algorithm progress
 * @param thisSpectrum :: [input] Spectrum number currently being fitted
 * @param totalSpectra :: [input] Total number of spectra to fit
 */
void MuonPhasesFromDetectors::reportProgress(const int thisSpectrum,
                                           const int totalSpectra) {
  double proportionDone = (double)thisSpectrum / (double)totalSpectra;
  std::ostringstream progMessage;
  progMessage << "Detector " << thisSpectrum + 1 << " of " << totalSpectra;
  this->progress(proportionDone, progMessage.str());
}

} // namespace Algorithms
} // namespace Mantid
