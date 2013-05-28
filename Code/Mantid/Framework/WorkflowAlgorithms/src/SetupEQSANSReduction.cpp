/*WIKI* 
Create a PropertyManager object setting the reduction options for EQSANS.
The property manager object is then added to the PropertyManagerDataService.

See [http://www.mantidproject.org/Reduction_for_HFIR_SANS SANS Reduction] documentation for details.

*WIKI*/
//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidWorkflowAlgorithms/SetupEQSANSReduction.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/ListValidator.h"
#include "MantidKernel/RebinParamsValidator.h"
#include "MantidKernel/EnabledWhenProperty.h"
#include "MantidKernel/VisibleWhenProperty.h"
#include "MantidAPI/FileProperty.h"
#include "MantidKernel/ArrayProperty.h"
#include "MantidAPI/AlgorithmProperty.h"
#include "MantidAPI/PropertyManagerDataService.h"
#include "MantidKernel/PropertyManager.h"
#include "Poco/NumberFormatter.h"

namespace Mantid
{
namespace WorkflowAlgorithms
{

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(SetupEQSANSReduction)

/// Sets documentation strings for this algorithm
void SetupEQSANSReduction::initDocs()
{
  this->setWikiSummary("Set up EQSANS SANS reduction options.");
  this->setOptionalMessage("Set up EQSANS SANS reduction options.");
}

using namespace Kernel;
using namespace API;
using namespace Geometry;

void SetupEQSANSReduction::init()
{
  // Load options
  std::string load_grp = "Load Options";
  declareProperty("UseConfigTOFCuts", false, "If true, the edges of the TOF distribution will be cut according to the configuration file");
  declareProperty("LowTOFCut", 0.0, "TOF value below which events will not be loaded into the workspace at load-time");
  declareProperty("HighTOFCut", 0.0, "TOF value above which events will not be loaded into the workspace at load-time");
  declareProperty("WavelengthStep", 0.1, "Wavelength steps to be used when rebinning the data before performing the reduction");
  declareProperty("UseConfigMask", false, "If true, the masking information found in the configuration file will be used");
  declareProperty("UseConfig", true, "If true, the best configuration file found will be used");
  declareProperty("CorrectForFlightPath", false, "If true, the TOF will be modified for the true flight path from the sample to the detector pixel");

  declareProperty("SkipTOFCorrection", false, "IF true, the EQSANS TOF correction will be skipped");
  declareProperty("PreserveEvents", true, "If true, the output workspace will be an event workspace");

  declareProperty("SampleDetectorDistance", EMPTY_DBL(), "Sample to detector distance to use (overrides meta data), in mm");
  declareProperty("SampleDetectorDistanceOffset", EMPTY_DBL(), "Offset to the sample to detector distance (use only when using the distance found in the meta data), in mm");

  declareProperty("SolidAngleCorrection", true, "If true, the solide angle correction will be applied to the data");
  declareProperty("DetectorTubes", false, "If true, the solid angle correction for tube detectors will be applied");

  // -- Define group --
  setPropertyGroup("UseConfigTOFCuts", load_grp);
  setPropertyGroup("LowTOFCut", load_grp);
  setPropertyGroup("HighTOFCut", load_grp);

  setPropertyGroup("WavelengthStep", load_grp);
  setPropertyGroup("UseConfigMask", load_grp);
  setPropertyGroup("UseConfig", load_grp);
  setPropertyGroup("CorrectForFlightPath", load_grp);

  setPropertyGroup("SkipTOFCorrection", load_grp);
  setPropertyGroup("PreserveEvents", load_grp);

  setPropertyGroup("SampleDetectorDistance", load_grp);
  setPropertyGroup("SampleDetectorDistanceOffset", load_grp);

  setPropertyGroup("SolidAngleCorrection", load_grp);
  setPropertyGroup("DetectorTubes", load_grp);

  // Beam center
  std::string center_grp = "Beam Center";
  std::vector<std::string> centerOptions;
  centerOptions.push_back("None");
  centerOptions.push_back("Value");
  centerOptions.push_back("DirectBeam");
  centerOptions.push_back("Scattering");

  declareProperty("BeamCenterMethod", "None",
      boost::make_shared<StringListValidator>(centerOptions),
      "Method for determining the data beam center");

  //declareProperty("FindBeamCenter", false, "If True, the beam center will be calculated");
  declareProperty("UseConfigBeam", false, "If True, the beam center will be taken from the config file");

  //    Option 1: Set beam center by hand
  declareProperty("BeamCenterX", EMPTY_DBL(), "Position of the beam center, in pixel");
  declareProperty("BeamCenterY", EMPTY_DBL(), "Position of the beam center, in pixel");
  setPropertySettings("BeamCenterX",
            new VisibleWhenProperty("BeamCenterMethod", IS_EQUAL_TO, "Value"));
  setPropertySettings("BeamCenterY",
            new VisibleWhenProperty("BeamCenterMethod", IS_EQUAL_TO, "Value"));

  //    Option 2: Find it (expose properties from FindCenterOfMass)
  declareProperty(new API::FileProperty("BeamCenterFile", "", API::FileProperty::OptionalLoad, "_event.nxs"),
      "The name of the input event Nexus file to load");
  setPropertySettings("BeamCenterFile",
            new VisibleWhenProperty("BeamCenterMethod", IS_NOT_EQUAL_TO, "None"));

  //declareProperty("Tolerance", EMPTY_DBL(), "Tolerance on the center of mass position between each iteration [m]. Default: 0.00125");
  auto positiveDouble = boost::make_shared<BoundedValidator<double> >();
  positiveDouble->setLower(0);
  //declareProperty("UseDirectBeamMethod", true, "If true, the direct beam method will be used");
  declareProperty("BeamRadius", EMPTY_DBL(),
      "Radius of the beam area used the exclude the beam when calculating "
      "the center of mass of the scattering pattern [pixels]. Default=3.0");
  setPropertySettings("BeamRadius",
            new VisibleWhenProperty("BeamCenterMethod", IS_EQUAL_TO, "Scattering"));

  // -- Define group --
  setPropertyGroup("BeamCenterMethod", center_grp);
  setPropertyGroup("UseConfigBeam", center_grp);
  setPropertyGroup("BeamCenterX", center_grp);
  setPropertyGroup("BeamCenterY", center_grp);
  setPropertyGroup("BeamCenterFile", center_grp);
  //setPropertyGroup("Tolerance", center_grp);
  //setPropertyGroup("UseDirectBeamMethod", center_grp);
  setPropertyGroup("BeamRadius", center_grp);


  // Normalisation
  std::string norm_grp = "Normalisation";
  std::vector<std::string> incidentBeamNormOptions;
  incidentBeamNormOptions.push_back("None");
  // The data will be normalised to the monitor counts
  incidentBeamNormOptions.push_back("Monitor");
  // The data will be normalised to the total charge and divided by the beam profile
  incidentBeamNormOptions.push_back("BeamProfileAndCharge");
  // The data will be normalised to the total charge only (no beam profile)
  incidentBeamNormOptions.push_back("Charge");
  this->declareProperty("Normalisation", "BeamProfileAndCharge",
      boost::make_shared<StringListValidator>(incidentBeamNormOptions),
      "Options for data normalisation");

  declareProperty("LoadMonitors", false, "If true, the monitor workspace will be loaded");
  //declareProperty("NormaliseToBeam", true, "If true, the data will be normalised to the total charge and divided by the beam profile");
  //declareProperty("NormaliseToMonitor", false, "If true, the data will be normalised to the monitor, otherwise the total charge will be used");
  declareProperty(new API::FileProperty("MonitorReferenceFile", "", API::FileProperty::OptionalLoad, "_event.nxs"),
      "The name of the beam monitor reference file used for normalisation");

  setPropertyGroup("Normalisation", norm_grp);
  setPropertyGroup("LoadMonitors", norm_grp);
  setPropertyGroup("MonitorReferenceFile", norm_grp);

  // Dark current
  declareProperty(new API::FileProperty("DarkCurrentFile", "", API::FileProperty::OptionalLoad, "_event.nxs"),
      "The name of the input event Nexus file to load as dark current.");

  // Sensitivity
  std::string eff_grp = "Sensitivity";
  declareProperty(new API::FileProperty("SensitivityFile", "", API::FileProperty::OptionalLoad, "_event.nxs"),
      "Flood field or sensitivity file.");
  declareProperty("MinEfficiency", EMPTY_DBL(), positiveDouble,
      "Minimum efficiency for a pixel to be considered (default: no minimum).");
  declareProperty("MaxEfficiency", EMPTY_DBL(), positiveDouble,
      "Maximum efficiency for a pixel to be considered (default: no maximum).");
  declareProperty("UseDefaultDC", true, "If true, the dark current subtracted from the sample data will also be subtracted from the flood field.");
  declareProperty(new API::FileProperty("SensitivityDarkCurrentFile", "", API::FileProperty::OptionalLoad, "_event.nxs"),
      "The name of the input file to load as dark current.");
  // - sensitivity beam center
  declareProperty("SensitivityBeamCenterMethod", "None",
      boost::make_shared<StringListValidator>(centerOptions),
      "Method for determining the sensitivity data beam center");

  //    Option 1: Set beam center by hand
  declareProperty("SensitivityBeamCenterX", EMPTY_DBL(),
      "Sensitivity beam center location in X [pixels]");
  setPropertySettings("SensitivityBeamCenterX",
            new VisibleWhenProperty("SensitivityBeamCenterMethod", IS_EQUAL_TO, "Value"));

  declareProperty("SensitivityBeamCenterY", EMPTY_DBL(),
      "Sensitivity beam center location in Y [pixels]");
  setPropertySettings("SensitivityBeamCenterY",
            new VisibleWhenProperty("SensitivityBeamCenterMethod", IS_EQUAL_TO, "Value"));

  //    Option 2: Find it (expose properties from FindCenterOfMass)
  declareProperty(new API::FileProperty("SensitivityBeamCenterFile", "",
      API::FileProperty::OptionalLoad, ".xml"),
      "The name of the input data file to load");
  setPropertySettings("SensitivityBeamCenterFile",
            new VisibleWhenProperty("SensitivityBeamCenterMethod", IS_NOT_EQUAL_TO, "None"));

  declareProperty("SensitivityBeamCenterRadius", EMPTY_DBL(),
      "Radius of the beam area used the exclude the beam when calculating "
      "the center of mass of the scattering pattern [pixels]. Default=3.0");
  setPropertySettings("SensitivityBeamCenterRadius",
            new VisibleWhenProperty("BeamCenterMethod", IS_EQUAL_TO, "Scattering"));

  declareProperty("OutputSensitivityWorkspace", "",
      "Name to give the sensitivity workspace");

  // -- Define group --
  setPropertyGroup("SensitivityFile", eff_grp);
  setPropertyGroup("MinEfficiency", eff_grp);
  setPropertyGroup("MaxEfficiency", eff_grp);
  setPropertyGroup("UseDefaultDC", eff_grp);
  setPropertyGroup("SensitivityDarkCurrentFile", eff_grp);
  setPropertyGroup("SensitivityBeamCenterMethod", eff_grp);
  setPropertyGroup("SensitivityBeamCenterX", eff_grp);
  setPropertyGroup("SensitivityBeamCenterY", eff_grp);
  setPropertyGroup("SensitivityBeamCenterFile", eff_grp);
  setPropertyGroup("SensitivityBeamCenterRadius", eff_grp);
  setPropertyGroup("OutputSensitivityWorkspace", eff_grp);

  // Transmission
  std::string trans_grp = "Transmission";
  std::vector<std::string> transOptions;
  transOptions.push_back("Value");
  transOptions.push_back("DirectBeam");
  declareProperty("TransmissionMethod", "Value",
      boost::make_shared<StringListValidator>(transOptions),
      "Transmission determination method");

  // - Transmission value entered by hand
  declareProperty("TransmissionValue", EMPTY_DBL(), positiveDouble,
      "Transmission value.");
  setPropertySettings("TransmissionValue",
            new VisibleWhenProperty("TransmissionMethod", IS_EQUAL_TO, "Value"));
  declareProperty("TransmissionError", EMPTY_DBL(), positiveDouble,
      "Transmission error.");
  setPropertySettings("TransmissionError",
            new VisibleWhenProperty("TransmissionMethod", IS_EQUAL_TO, "Value"));

  // - Direct beam method transmission calculation
  declareProperty("TransmissionBeamRadius", 3.0,
      "Radius of the beam area used to compute the transmission [pixels]");
  setPropertySettings("TransmissionBeamRadius",
            new VisibleWhenProperty("TransmissionMethod", IS_EQUAL_TO, "DirectBeam"));
  declareProperty(new API::FileProperty("TransmissionSampleDataFile", "",
      API::FileProperty::OptionalLoad, ".xml"),
      "Sample data file for transmission calculation");
  setPropertySettings("TransmissionSampleDataFile",
            new VisibleWhenProperty("TransmissionMethod", IS_EQUAL_TO, "DirectBeam"));
  declareProperty(new API::FileProperty("TransmissionEmptyDataFile", "",
      API::FileProperty::OptionalLoad, ".xml"),
      "Empty data file for transmission calculation");
  setPropertySettings("TransmissionEmptyDataFile",
            new VisibleWhenProperty("TransmissionMethod", IS_EQUAL_TO, "DirectBeam"));
  declareProperty("FitFramesTogether", false,
      "If true, the two frames will be fit together");
  setPropertySettings("FitFramesTogether",
            new VisibleWhenProperty("TransmissionMethod", IS_EQUAL_TO, "DirectBeam"));

  // - transmission beam center
  declareProperty("TransmissionBeamCenterMethod", "None",
      boost::make_shared<StringListValidator>(centerOptions),
      "Method for determining the transmission data beam center");
  setPropertySettings("TransmissionBeamCenterMethod",
            new VisibleWhenProperty("TransmissionMethod", IS_EQUAL_TO, "DirectBeam"));

  //    Option 1: Set beam center by hand
  declareProperty("TransmissionBeamCenterX", EMPTY_DBL(),
      "Transmission beam center location in X [pixels]");
  setPropertySettings("TransmissionBeamCenterX",
            new VisibleWhenProperty("TransmissionMethod", IS_EQUAL_TO, "DirectBeam"));
  declareProperty("TransmissionBeamCenterY", EMPTY_DBL(),
      "Transmission beam center location in Y [pixels]");
  setPropertySettings("TransmissionBeamCenterY",
            new VisibleWhenProperty("TransmissionMethod", IS_EQUAL_TO, "DirectBeam"));

  //    Option 2: Find it (expose properties from FindCenterOfMass)
  declareProperty(new API::FileProperty("TransmissionBeamCenterFile", "",
      API::FileProperty::OptionalLoad, ".xml"),
      "The name of the input data file to load");
  setPropertySettings("TransmissionBeamCenterFile",
            new VisibleWhenProperty("TransmissionMethod", IS_EQUAL_TO, "DirectBeam"));

  declareProperty(new API::FileProperty("TransmissionDarkCurrentFile", "", API::FileProperty::OptionalLoad, ".xml"),
      "The name of the input data file to load as transmission dark current.");
  setPropertySettings("TransmissionDarkCurrentFile",
            new VisibleWhenProperty("TransmissionMethod", IS_NOT_EQUAL_TO, "Value"));

  declareProperty("TransmissionUseSampleDC", true,
      "If true, the sample dark current will be used IF a dark current file is"
      "not set.");
  setPropertySettings("TransmissionUseSampleDC",
            new VisibleWhenProperty("TransmissionMethod", IS_NOT_EQUAL_TO, "Value"));

  declareProperty("ThetaDependentTransmission", true,
      "If true, a theta-dependent transmission correction will be applied.");

  // -- Define group --
  setPropertyGroup("TransmissionMethod", trans_grp);
  setPropertyGroup("TransmissionValue", trans_grp);
  setPropertyGroup("TransmissionError", trans_grp);
  setPropertyGroup("TransmissionBeamRadius", trans_grp);
  setPropertyGroup("TransmissionSampleDataFile", trans_grp);
  setPropertyGroup("TransmissionEmptyDataFile", trans_grp);
  setPropertyGroup("FitFramesTogether", trans_grp);
  setPropertyGroup("TransmissionBeamCenterMethod", trans_grp);
  setPropertyGroup("TransmissionBeamCenterX", trans_grp);
  setPropertyGroup("TransmissionBeamCenterY", trans_grp);
  setPropertyGroup("TransmissionBeamCenterFile", trans_grp);

  setPropertyGroup("TransmissionDarkCurrentFile", trans_grp);
  setPropertyGroup("TransmissionUseSampleDC", trans_grp);
  setPropertyGroup("ThetaDependentTransmission", trans_grp);

  // Background options
  std::string bck_grp = "Background";
  declareProperty("BackgroundFiles", "", "Background data files");
  declareProperty("BckTransmissionMethod", "Value",
      boost::make_shared<StringListValidator>(transOptions),
      "Transmission determination method");

  // - Transmission value entered by hand
  declareProperty("BckTransmissionValue", EMPTY_DBL(), positiveDouble,
      "Transmission value.");
  setPropertySettings("BckTransmissionValue",
            new VisibleWhenProperty("BckTransmissionMethod", IS_EQUAL_TO, "Value"));

  declareProperty("BckTransmissionError", EMPTY_DBL(), positiveDouble,
      "Transmission error.");
  setPropertySettings("BckTransmissionError",
            new VisibleWhenProperty("BckTransmissionMethod", IS_EQUAL_TO, "Value"));

  // - Direct beam method transmission calculation
  declareProperty("BckTransmissionBeamRadius", 3.0,
      "Radius of the beam area used to compute the transmission [pixels]");
  setPropertySettings("BckTransmissionBeamRadius",
            new VisibleWhenProperty("BckTransmissionMethod", IS_EQUAL_TO, "DirectBeam"));
  declareProperty(new API::FileProperty("BckTransmissionSampleDataFile", "",
      API::FileProperty::OptionalLoad, ".xml"),
      "Sample data file for transmission calculation");
  setPropertySettings("BckTransmissionSampleDataFile",
            new VisibleWhenProperty("BckTransmissionMethod", IS_EQUAL_TO, "DirectBeam"));
  declareProperty(new API::FileProperty("BckTransmissionEmptyDataFile", "",
      API::FileProperty::OptionalLoad, ".xml"),
      "Empty data file for transmission calculation");
  setPropertySettings("BckTransmissionEmptyDataFile",
            new VisibleWhenProperty("BckTransmissionMethod", IS_EQUAL_TO, "DirectBeam"));
  declareProperty("BckFitFramesTogether", false,
      "If true, the two frames will be fit together");
  setPropertySettings("BckFitFramesTogether",
            new VisibleWhenProperty("BckTransmissionMethod", IS_EQUAL_TO, "DirectBeam"));

  // - transmission beam center
  declareProperty("BckTransmissionBeamCenterMethod", "None",
      boost::make_shared<StringListValidator>(centerOptions),
      "Method for determining the transmission data beam center");
  setPropertySettings("BckTransmissionBeamCenterMethod",
            new VisibleWhenProperty("BckTransmissionMethod", IS_EQUAL_TO, "DirectBeam"));
  //    Option 1: Set beam center by hand
  declareProperty("BckTransmissionBeamCenterX", EMPTY_DBL(),
      "Transmission beam center location in X [pixels]");
  setPropertySettings("BckTransmissionBeamCenterX",
            new VisibleWhenProperty("BckTransmissionMethod", IS_EQUAL_TO, "DirectBeam"));
  declareProperty("BckTransmissionBeamCenterY", EMPTY_DBL(),
      "Transmission beam center location in Y [pixels]");
  //    Option 2: Find it (expose properties from FindCenterOfMass)
  setPropertySettings("BckTransmissionBeamCenterY",
            new VisibleWhenProperty("BckTransmissionMethod", IS_EQUAL_TO, "DirectBeam"));
  declareProperty(new API::FileProperty("BckTransmissionBeamCenterFile", "",
      API::FileProperty::OptionalLoad, ".xml"),
      "The name of the input data file to load");
  setPropertySettings("BckTransmissionBeamCenterFile",
            new VisibleWhenProperty("BckTransmissionMethod", IS_EQUAL_TO, "DirectBeam"));

  declareProperty(new API::FileProperty("BckTransmissionDarkCurrentFile", "", API::FileProperty::OptionalLoad, ".xml"),
      "The name of the input data file to load as background transmission dark current.");
  setPropertySettings("BckTransmissionDarkCurrentFile",
            new VisibleWhenProperty("BckTransmissionMethod", IS_EQUAL_TO, "BeamSpreader"));

  declareProperty("BckThetaDependentTransmission", true,
      "If true, a theta-dependent transmission correction will be applied.");

  setPropertyGroup("BackgroundFiles", bck_grp);
  setPropertyGroup("BckTransmissionMethod", bck_grp);
  setPropertyGroup("BckTransmissionValue", bck_grp);
  setPropertyGroup("BckTransmissionError", bck_grp);
  setPropertyGroup("BckTransmissionBeamRadius", bck_grp);
  setPropertyGroup("BckTransmissionSampleDataFile", bck_grp);
  setPropertyGroup("BckTransmissionEmptyDataFile", bck_grp);
  setPropertyGroup("BckTransmissionBeamCenterMethod", bck_grp);
  setPropertyGroup("BckTransmissionBeamCenterX", bck_grp);
  setPropertyGroup("BckTransmissionBeamCenterY", bck_grp);
  setPropertyGroup("BckTransmissionBeamCenterFile", bck_grp);
  setPropertyGroup("BckTransmissionDarkCurrentFile", bck_grp);
  setPropertyGroup("BckThetaDependentTransmission", bck_grp);

  // Geometry correction
  declareProperty("SampleThickness", EMPTY_DBL(), "Sample thickness [cm]");

  // Masking
  std::string mask_grp = "Mask";
  declareProperty(new ArrayProperty<int>("MaskedDetectorList"),
      "List of detector IDs to be masked");
  declareProperty(new ArrayProperty<int>("MaskedEdges"),
      "Number of pixels to mask on the edges: X-low, X-high, Y-low, Y-high");
  std::vector<std::string> maskOptions;
  maskOptions.push_back("None");
  maskOptions.push_back("Front");
  maskOptions.push_back("Back");
  declareProperty("MaskedSide", "None",
      boost::make_shared<StringListValidator>(maskOptions),
      "Mask one side of the detector");

  setPropertyGroup("MaskedDetectorList", mask_grp);
  setPropertyGroup("MaskedEdges", mask_grp);
  setPropertyGroup("MaskedSide", mask_grp);

  // Absolute scale
  std::string abs_scale_grp = "Absolute Scale";
  std::vector<std::string> scaleOptions;
  scaleOptions.push_back("None");
  scaleOptions.push_back("Value");
  scaleOptions.push_back("ReferenceData");
  declareProperty("AbsoluteScaleMethod", "None",
      boost::make_shared<StringListValidator>(scaleOptions),
      "Absolute scale correction method");
  declareProperty("AbsoluteScalingFactor", 1.0, "Absolute scaling factor");
  setPropertySettings("AbsoluteScalingFactor",
      new VisibleWhenProperty("AbsoluteScaleMethod", IS_EQUAL_TO, "Value"));

  declareProperty(new API::FileProperty("AbsoluteScalingReferenceFilename", "",
      API::FileProperty::OptionalLoad, ".xml"));
  setPropertySettings("AbsoluteScalingReferenceFilename",
            new VisibleWhenProperty("AbsoluteScaleMethod", IS_EQUAL_TO, "ReferenceData"));
  declareProperty("AbsoluteScalingBeamDiameter", 0.0,
      "Beamstop diameter for computing the absolute scale factor [mm]. "
      "Read from file if not supplied.");
  setPropertySettings("AbsoluteScalingBeamDiameter",
      new VisibleWhenProperty("AbsoluteScaleMethod", IS_EQUAL_TO, "ReferenceData"));
  declareProperty("AbsoluteScalingAttenuatorTrans", 1.0,
      "Attenuator transmission value for computing the absolute scale factor");
  setPropertySettings("AbsoluteScalingAttenuatorTrans",
      new VisibleWhenProperty("AbsoluteScaleMethod", IS_EQUAL_TO, "ReferenceData"));
  declareProperty("AbsoluteScalingApplySensitivity", false,
      "Apply sensitivity correction to the reference data "
      "when computing the absolute scale factor");
  setPropertySettings("AbsoluteScalingApplySensitivity",
      new VisibleWhenProperty("AbsoluteScaleMethod", IS_EQUAL_TO, "ReferenceData"));

  setPropertyGroup("AbsoluteScaleMethod", abs_scale_grp);
  setPropertyGroup("AbsoluteScalingFactor", abs_scale_grp);
  setPropertyGroup("AbsoluteScalingReferenceFilename", abs_scale_grp);
  setPropertyGroup("AbsoluteScalingBeamDiameter", abs_scale_grp);
  setPropertyGroup("AbsoluteScalingAttenuatorTrans", abs_scale_grp);
  setPropertyGroup("AbsoluteScalingApplySensitivity", abs_scale_grp);


  // Setup the version 1 reducer
  declareProperty("SetupReducer",false, "If true, a Reducer object will be created");

  // I(Q) calculation
  std::string iq1d_grp = "I(q) Calculation";
  declareProperty("DoAzimuthalAverage", true);
  auto positiveInt = boost::make_shared<BoundedValidator<int> >();
  positiveInt->setLower(0);
  declareProperty("IQNumberOfBins", 100, positiveInt,
                  "Number of I(q) bins when binning is not specified");
  declareProperty("IQLogBinning", false,
                  "I(q) log binning when binning is not specified");
  declareProperty("IQIndependentBinning", true,
                  "If true and frame skipping is used, each frame will have its own binning");
  declareProperty("IQScaleResults", true,
                  "If true and frame skipping is used, frame 1 will be scaled to frame 2");
  declareProperty("ComputeResolution", false,
                  "If true the Q resolution will be computed");
  declareProperty("SampleApertureDiameter", 10.0,
                  "Sample aperture diameter [mm]");

  declareProperty("Do2DReduction", true);
  declareProperty("IQ2DNumberOfBins", 100, positiveInt,
                  "Number of I(qx,qy) bins.");

  // -- Define group --
  setPropertyGroup("DoAzimuthalAverage", iq1d_grp);
  setPropertyGroup("IQNumberOfBins", iq1d_grp);
  setPropertyGroup("IQLogBinning", iq1d_grp);
  setPropertyGroup("IQIndependentBinning", iq1d_grp);
  setPropertyGroup("IQScaleResults", iq1d_grp);
  setPropertyGroup("ComputeResolution", iq1d_grp);
  setPropertyGroup("SampleApertureDiameter", iq1d_grp);
  setPropertyGroup("Do2DReduction", iq1d_grp);
  setPropertyGroup("IQ2DNumberOfBins", iq1d_grp);

  // Outputs
  declareProperty("ProcessInfo","", "Additional process information");
  declareProperty("OutputDirectory", "", "Directory to put the output files in");
  declareProperty("OutputMessage","",Direction::Output);
  declareProperty("ReductionProperties","__sans_reduction_properties", Direction::Input);
}

void SetupEQSANSReduction::exec()
{
  // Reduction property manager
  const std::string reductionManagerName = getProperty("ReductionProperties");
  if (reductionManagerName.size()==0)
  {
    g_log.error() << "ERROR: Reduction Property Manager name is empty" << std::endl;
    return;
  }
  boost::shared_ptr<PropertyManager> reductionManager = boost::make_shared<PropertyManager>();
  PropertyManagerDataService::Instance().addOrReplace(reductionManagerName, reductionManager);

  // Store name of the instrument
  reductionManager->declareProperty(new PropertyWithValue<std::string>("InstrumentName", "EQSANS") );

  // Store additional (and optional) process information
  const std::string processInfo = getProperty("ProcessInfo");
  reductionManager->declareProperty(new PropertyWithValue<std::string>("ProcessInfo", processInfo));

  // Store the output directory
  const std::string outputDirectory = getProperty("OutputDirectory");
  reductionManager->declareProperty(new PropertyWithValue<std::string>("OutputDirectory", outputDirectory));

  // Store normalization algorithm
  const std::string normalization = getProperty("Normalisation");
  bool loadMonitors = getProperty("LoadMonitors");
  const std::string monitorRefFile = getPropertyValue("MonitorReferenceFile");
  // If we normalize to monitor, force the loading of monitor data
  IAlgorithm_sptr normAlg = createChildAlgorithm("EQSANSNormalise");

  if (boost::contains(normalization, "BeamProfileAndCharge"))
  {
    normAlg->setProperty("NormaliseToBeam", true);
    normAlg->setProperty("BeamSpectrumFile", monitorRefFile);
  }
  else if (boost::contains(normalization, "Charge"))
  {
    normAlg->setProperty("NormaliseToBeam", false);
  }
  else if (boost::contains(normalization, "Monitor"))
  {
    loadMonitors = true;
    if (monitorRefFile.size()==0)
    {
      g_log.error() << "ERROR: normalize-to-monitor was turned ON but no reference data was selected" << std::endl;
    }
    normAlg->setProperty("NormaliseToMonitor", true);
    normAlg->setProperty("BeamSpectrumFile", monitorRefFile);
  }
  normAlg->setPropertyValue("ReductionProperties", reductionManagerName);
  AlgorithmProperty *algProp = new AlgorithmProperty("NormaliseAlgorithm");
  algProp->setValue(normAlg->toString());
  reductionManager->declareProperty(algProp);

  // Load algorithm
  IAlgorithm_sptr loadAlg = createChildAlgorithm("EQSANSLoad");
  const bool useConfigBeam = getProperty("UseConfigBeam");
  loadAlg->setProperty("UseConfigBeam", useConfigBeam);
  const bool useConfigTOFCuts = getProperty("UseConfigTOFCuts");
  loadAlg->setProperty("UseConfigTOFCuts", useConfigTOFCuts);
  if (!useConfigTOFCuts)
  {
    const double lowTOFCut = getProperty("LowTOFCut");
    const double highTOFCut = getProperty("HighTOFCut");
    loadAlg->setProperty("LowTOFCut", lowTOFCut);
    loadAlg->setProperty("HighTOFCut", highTOFCut);
  }

  const bool skipTOFCorrection = getProperty("SkipTOFCorrection");
  loadAlg->setProperty("SkipTOFCorrection", skipTOFCorrection);

  const bool correctForFlightPath = getProperty("CorrectForFlightPath");
  loadAlg->setProperty("CorrectForFlightPath", correctForFlightPath);

  const bool preserveEvents = getProperty("PreserveEvents");
  loadAlg->setProperty("PreserveEvents", preserveEvents);
  loadAlg->setProperty("LoadMonitors", loadMonitors);

  const double sdd = getProperty("SampleDetectorDistance");
  loadAlg->setProperty("SampleDetectorDistance", sdd);
  const double sddOffset = getProperty("SampleDetectorDistanceOffset");
  loadAlg->setProperty("SampleDetectorDistanceOffset", sddOffset);
  const double wlStep = getProperty("WavelengthStep");
  loadAlg->setProperty("WavelengthStep", wlStep);

  const bool useConfig = getProperty("UseConfig");
  loadAlg->setProperty("UseConfig", useConfig);
  const bool useConfigMask = getProperty("UseConfigMask");
  loadAlg->setProperty("UseConfigMask", useConfigMask);
  algProp = new AlgorithmProperty("LoadAlgorithm");
  algProp->setValue(loadAlg->toString());
  reductionManager->declareProperty(algProp);

  // Store dark current algorithm
  const std::string darkCurrentFile = getPropertyValue("DarkCurrentFile");
  if (darkCurrentFile.size() > 0)
  {
    IAlgorithm_sptr darkAlg = createChildAlgorithm("EQSANSDarkCurrentSubtraction");
    darkAlg->setProperty("Filename", darkCurrentFile);
    darkAlg->setProperty("OutputDarkCurrentWorkspace", "");
    darkAlg->setPropertyValue("ReductionProperties", reductionManagerName);
    algProp = new AlgorithmProperty("DarkCurrentAlgorithm");
    algProp->setValue(darkAlg->toString());
    reductionManager->declareProperty(algProp);
  }

  // Store default dark current algorithm
  IAlgorithm_sptr darkDefaultAlg = createChildAlgorithm("EQSANSDarkCurrentSubtraction");
  darkDefaultAlg->setProperty("OutputDarkCurrentWorkspace", "");
  darkDefaultAlg->setPropertyValue("ReductionProperties", reductionManagerName);
  algProp = new AlgorithmProperty("DefaultDarkCurrentAlgorithm");
  algProp->setValue(darkDefaultAlg->toString());
  reductionManager->declareProperty(algProp);

  // Solid angle correction
  const bool solidAngleCorrection = getProperty("SolidAngleCorrection");
  if (solidAngleCorrection)
  {
    const bool detectorTubes = getProperty("DetectorTubes");
    IAlgorithm_sptr solidAlg = createChildAlgorithm("SANSSolidAngleCorrection");
    solidAlg->setProperty("DetectorTubes", detectorTubes);
    algProp = new AlgorithmProperty("SANSSolidAngleCorrection");
    algProp->setValue(solidAlg->toString());
    reductionManager->declareProperty(algProp);
  }

  // Beam center
  const double beamCenterX = getProperty("BeamCenterX");
  const double beamCenterY = getProperty("BeamCenterY");
  const std::string centerMethod = getPropertyValue("BeamCenterMethod");

  // Beam center option for transmission data
  if (boost::iequals(centerMethod, "Value"))
  {
    if(!isEmpty(beamCenterX) && !isEmpty(beamCenterY))
    {
      reductionManager->declareProperty(new PropertyWithValue<double>("LatestBeamCenterX", beamCenterX) );
      reductionManager->declareProperty(new PropertyWithValue<double>("LatestBeamCenterY", beamCenterY) );
    }
  }
  else if (!boost::iequals(centerMethod, "None"))
  {
    bool useDirectBeamMethod = true;
    if (!boost::iequals(centerMethod, "DirectBeam")) useDirectBeamMethod = false;
    const std::string beamCenterFile = getProperty("BeamCenterFile");
    if (beamCenterFile.size()>0)
    {
      const double beamRadius = getProperty("BeamRadius");

      IAlgorithm_sptr ctrAlg = createChildAlgorithm("SANSBeamFinder");
      ctrAlg->setProperty("Filename", beamCenterFile);
      ctrAlg->setProperty("UseDirectBeamMethod", useDirectBeamMethod);
      if (!isEmpty(beamRadius)) ctrAlg->setProperty("BeamRadius", beamRadius);
      ctrAlg->setPropertyValue("ReductionProperties", reductionManagerName);

      AlgorithmProperty *algProp = new AlgorithmProperty("SANSBeamFinderAlgorithm");
      algProp->setValue(ctrAlg->toString());
      reductionManager->declareProperty(algProp);
    } else {
      g_log.error() << "ERROR: Beam center determination was required"
          " but no file was provided" << std::endl;
    }
  }

  // Sensitivity correction, transmission and background
  setupSensitivity(reductionManager);
  setupTransmission(reductionManager);
  setupBackground(reductionManager);

  // Geometry correction
  const double thickness = getProperty("SampleThickness");
  if (!isEmpty(thickness))
  {
    IAlgorithm_sptr thickAlg = createChildAlgorithm("NormaliseByThickness");
    thickAlg->setProperty("SampleThickness", thickness);

    algProp = new AlgorithmProperty("GeometryAlgorithm");
    algProp->setValue(thickAlg->toString());
    reductionManager->declareProperty(algProp);
  }

  // Mask
  const std::string maskDetList = getPropertyValue("MaskedDetectorList");
  const std::string maskEdges = getPropertyValue("MaskedEdges");
  const std::string maskSide = getProperty("MaskedSide");

  IAlgorithm_sptr maskAlg = createChildAlgorithm("SANSMask");
  // The following is broken, try PropertyValue
  maskAlg->setPropertyValue("Facility", "SNS");
  maskAlg->setPropertyValue("MaskedDetectorList", maskDetList);
  maskAlg->setPropertyValue("MaskedEdges", maskEdges);
  maskAlg->setProperty("MaskedSide", maskSide);
  algProp = new AlgorithmProperty("MaskAlgorithm");
  algProp->setValue(maskAlg->toString());
  reductionManager->declareProperty(algProp);

  // Absolute scaling
  const std::string absScaleMethod = getProperty("AbsoluteScaleMethod");
  if (boost::iequals(absScaleMethod, "Value"))
  {
    const double absScaleFactor = getProperty("AbsoluteScalingFactor");

    IAlgorithm_sptr absAlg = createChildAlgorithm("SANSAbsoluteScale");
    absAlg->setProperty("Method", absScaleMethod);
    absAlg->setProperty("ScalingFactor", absScaleFactor);
    absAlg->setPropertyValue("ReductionProperties", reductionManagerName);
    algProp = new AlgorithmProperty("AbsoluteScaleAlgorithm");
    algProp->setValue(absAlg->toString());
    reductionManager->declareProperty(algProp);
  }
  else if (boost::iequals(absScaleMethod, "ReferenceData"))
  {
    const std::string absRefFile = getPropertyValue("AbsoluteScalingReferenceFilename");
    const double beamDiam = getProperty("AbsoluteScalingBeamDiameter");
    const double attTrans = getProperty("AbsoluteScalingAttenuatorTrans");
    const bool applySensitivity = getProperty("AbsoluteScalingApplySensitivity");

    IAlgorithm_sptr absAlg = createChildAlgorithm("SANSAbsoluteScale");
    absAlg->setProperty("Method", absScaleMethod);
    absAlg->setProperty("ReferenceDataFilename", absRefFile);
    absAlg->setProperty("BeamstopDiameter", beamDiam);
    absAlg->setProperty("AttenuatorTransmission", attTrans);
    absAlg->setProperty("ApplySensitivity", applySensitivity);
    absAlg->setPropertyValue("ReductionProperties", reductionManagerName);
    algProp = new AlgorithmProperty("AbsoluteScaleAlgorithm");
    algProp->setValue(absAlg->toString());
    reductionManager->declareProperty(algProp);
  }

  // Azimuthal averaging
  const bool doAveraging = getProperty("DoAzimuthalAverage");
  if (doAveraging)
  {
    const std::string nBins = getPropertyValue("IQNumberOfBins");
    const bool logBinning = getProperty("IQLogBinning");
    const double sampleApert = getProperty("SampleApertureDiameter");
    const bool computeResolution = getProperty("ComputeResolution");
    const bool indepBinning = getProperty("IQIndependentBinning");
    const bool scaleResults = getProperty("IQScaleResults");

    IAlgorithm_sptr iqAlg = createChildAlgorithm("EQSANSAzimuthalAverage1D");
    iqAlg->setPropertyValue("NumberOfBins", nBins);
    iqAlg->setProperty("LogBinning", logBinning);
    iqAlg->setProperty("ScaleResults", scaleResults);
    iqAlg->setProperty("ComputeResolution", computeResolution);
    iqAlg->setProperty("IndependentBinning", indepBinning);
    iqAlg->setProperty("SampleApertureDiameter", sampleApert);
    iqAlg->setPropertyValue("ReductionProperties", reductionManagerName);

    algProp = new AlgorithmProperty("IQAlgorithm");
    algProp->setValue(iqAlg->toString());
    reductionManager->declareProperty(algProp);
  }

  // 2D reduction
  const bool do2DReduction = getProperty("Do2DReduction");
  if (do2DReduction)
  {
    const std::string n_bins = getPropertyValue("IQ2DNumberOfBins");
    IAlgorithm_sptr iqAlg = createChildAlgorithm("EQSANSQ2D");
    iqAlg->setPropertyValue("NumberOfBins", n_bins);
    algProp = new AlgorithmProperty("IQXYAlgorithm");
    algProp->setValue(iqAlg->toString());
    reductionManager->declareProperty(algProp);
  }
  setPropertyValue("OutputMessage", "EQSANS reduction options set");

  // Create a python reduction singleton as needed
  const bool setupReducer = getProperty("SetupReducer");
  if (setupReducer) initializeReduction(reductionManager);
}

/*
 * For backward compatibility, we have the option of creating a
 * python ReductionSingleton object.
 */
void SetupEQSANSReduction::initializeReduction(boost::shared_ptr<PropertyManager> reductionManager)
{
  // Write the Reducer python script to be executed
  std::string script = "import reduction.instruments.sans.sns_command_interface as cmd\n";
  const std::string reductionManagerName = getProperty("ReductionProperties");
  const bool preserveEvents = getProperty("PreserveEvents");

  //  - beam center
  double center_x = 0.0;
  double center_y = 0.0;
  if (reductionManager->existsProperty("LatestBeamCenterX")
      && reductionManager->existsProperty("LatestBeamCenterY"))
  {
    center_x = reductionManager->getProperty("LatestBeamCenterX");
    center_y = reductionManager->getProperty("LatestBeamCenterY");
  }
  else
    throw std::runtime_error("EQSANSReduce not yet compatible with beam finder: enter beam center coordinates");

  if (preserveEvents)
    script += "cmd.EQSANS(True, \"" + reductionManagerName + "\")\n";
  else
    script += "cmd.EQSANS(False, \"" + reductionManagerName + "\")\n";
  script += "cmd.SetBeamCenter(" + Poco::NumberFormatter::format(center_x, 2)
    + ", " + Poco::NumberFormatter::format(center_y, 2) + ")\n";

  //  - sensitivity file
  if (reductionManager->existsProperty("SensitivityAlgorithm"))
  {
    IAlgorithm_sptr effAlg = reductionManager->getProperty("SensitivityAlgorithm");
    const std::string fileName = effAlg->getPropertyValue("Filename");
    if (fileName.size()>0)
      script += "cmd.SensitivityCorrection(\"" + fileName + "\")\n";
  }

  //  - load options
  const bool useConfig = getProperty("UseConfig");
  if (useConfig)
    script += "cmd.UseConfig(True)\n";
  else
    script += "cmd.UseConfig(False)\n";

  if (reductionManager->existsProperty("LoadAlgorithm"))
  {
    IAlgorithm_sptr loadAlg = reductionManager->getProperty("LoadAlgorithm");

    // Correct for flight path?
    const bool tofCorr = loadAlg->getProperty("CorrectForFlightPath");
    if (tofCorr)
      script += "cmd.PerformFlightPathCorrection(True)\n";
    else
      script += "cmd.PerformFlightPathCorrection(False)\n";

    // Use TOF cut from config file?
    const bool confTOF = loadAlg->getProperty("UseConfigTOFCuts");
    if (confTOF)
      script += "cmd.UseConfigTOFTailsCutoff(use_config=True)\n";
    else
    {
      script += "cmd.UseConfigTOFTailsCutoff(use_config=False)\n";
      // Manual TOF cut
      const double lowTOF = loadAlg->getProperty("LowTOFCut");
      const double highTOF = loadAlg->getProperty("HighTOFCut");
      script += "cmd.SetTOFTailsCutoff(low_cut=" + Poco::NumberFormatter::format(lowTOF, 2)
          + ", high_cut=" + Poco::NumberFormatter::format(highTOF, 2) + ")\n";
    }

    // Use config mask?
    const bool confMask = loadAlg->getProperty("UseConfigMask");
    if (confMask)
      script += "cmd.UseConfigMask(use_config=True)\n";
    else
      script += "cmd.UseConfigMask(use_config=False)\n";
  }

  //  - Solid angle correction
  const bool solidAngle = getProperty("SolidAngleCorrection");
  if (solidAngle)
    script += "cmd.SolidAngle()\n";
  else
    script += "cmd.NoSolidAngle()\n";

  //  - Dark current
  const std::string darkCurrentFile = getPropertyValue("DarkCurrentFile");
  if (darkCurrentFile.size() > 0)
  {
    script += "cmd.DarkCurrent(\"" + darkCurrentFile + "\")\n";
  }

  //  - Normalization options
  if (reductionManager->existsProperty("NormaliseAlgorithm"))
  {
    IAlgorithm_sptr normAlg = reductionManager->getProperty("NormaliseAlgorithm");
    const bool normaliseToBeam = normAlg->getProperty("NormaliseToBeam");
    const bool normaliseToMonitor = normAlg->getProperty("NormaliseToMonitor");
    const std::string fileName = normAlg->getPropertyValue("BeamSpectrumFile");
    if (normaliseToMonitor)
    {
      script += "cmd.BeamMonitorNormalization(\"" + fileName + "\")\n";
    } else {
      if (normaliseToBeam)
        script += "cmd.TotalChargeNormalization(normalize_to_beam=True)\n";
      else
        script += "cmd.TotalChargeNormalization(normalize_to_beam=False)\n";
    }
  }

  //  - Transmission
  const double trans = getProperty("TransmissionValue");
  const bool thetaDependent = getProperty("ThetaDependentTransmission");
  if (isEmpty(trans))
  {
    const std::string directBeam = getPropertyValue("TransmissionDirectBeam");
    const std::string emptyBeam = getPropertyValue("TransmissionEmptyBeam");

    script += "cmd.DirectBeamTransmission(\"" + directBeam
        + "\", \"" + emptyBeam + "\")\n";
    if (thetaDependent)
      script += "cmd.ThetaDependentTransmission(True)\n";
    else
      script += "cmd.ThetaDependentTransmission(False)\n";
  }
  else
  {
    if (thetaDependent)
      script += "cmd.SetTransmission(" + Poco::NumberFormatter::format(trans, 2) + ", 0.0, True)\n";
    else
      script += "cmd.SetTransmission(" + Poco::NumberFormatter::format(trans, 2) + ", 0.0, False)\n";
  }
  script += "cmd.ReductionSingleton().set_azimuthal_averager(None)\n";

  g_log.information() << "Reducer script:\n" << script << std::endl;

  // Run a snippet of python
  IAlgorithm_sptr alg = this->createChildAlgorithm("RunOldPythonScript");
  alg->setLogging(false);
  alg->setPropertyValue("Code", script);
  alg->execute();
}

void SetupEQSANSReduction::setupSensitivity(boost::shared_ptr<PropertyManager> reductionManager)
{
  const std::string reductionManagerName = getProperty("ReductionProperties");

  const std::string sensitivityFile = getPropertyValue("SensitivityFile");
  if (sensitivityFile.size() > 0)
  {
    const bool useSampleDC = getProperty("UseDefaultDC");
    const std::string sensitivityDarkCurrentFile = getPropertyValue("SensitivityDarkCurrentFile");
    const std::string outputSensitivityWS = getPropertyValue("OutputSensitivityWorkspace");
    const double minEff = getProperty("MinEfficiency");
    const double maxEff = getProperty("MaxEfficiency");
    const double sensitivityBeamCenterX = getProperty("SensitivityBeamCenterX");
    const double sensitivityBeamCenterY = getProperty("SensitivityBeamCenterY");

    IAlgorithm_sptr effAlg = createChildAlgorithm("SANSSensitivityCorrection");
    effAlg->setProperty("Filename", sensitivityFile);
    effAlg->setProperty("UseSampleDC", useSampleDC);
    effAlg->setProperty("DarkCurrentFile", sensitivityDarkCurrentFile);
    effAlg->setProperty("MinEfficiency", minEff);
    effAlg->setProperty("MaxEfficiency", maxEff);

    // Beam center option for sensitivity data
    const std::string centerMethod = getPropertyValue("SensitivityBeamCenterMethod");
    if (boost::iequals(centerMethod, "Value"))
    {
      if (!isEmpty(sensitivityBeamCenterX) &&
          !isEmpty(sensitivityBeamCenterY))
      {
        effAlg->setProperty("BeamCenterX", sensitivityBeamCenterX);
        effAlg->setProperty("BeamCenterY", sensitivityBeamCenterY);
      }
    }
    else if (boost::iequals(centerMethod, "DirectBeam") ||
        boost::iequals(centerMethod, "Scattering"))
    {
      const std::string beamCenterFile = getProperty("SensitivityBeamCenterFile");
      const double sensitivityBeamRadius = getProperty("SensitivityBeamCenterRadius");
      bool useDirectBeam = boost::iequals(centerMethod, "DirectBeam");
      if (beamCenterFile.size()>0)
       {
         IAlgorithm_sptr ctrAlg = createChildAlgorithm("SANSBeamFinder");
         ctrAlg->setProperty("Filename", beamCenterFile);
         ctrAlg->setProperty("UseDirectBeamMethod", useDirectBeam);
         ctrAlg->setProperty("PersistentCorrection", false);
         if (useDirectBeam && !isEmpty(sensitivityBeamRadius))
           ctrAlg->setProperty("BeamRadius", sensitivityBeamRadius);
         ctrAlg->setPropertyValue("ReductionProperties", reductionManagerName);

         AlgorithmProperty *algProp = new AlgorithmProperty("SensitivityBeamCenterAlgorithm");
         algProp->setValue(ctrAlg->toString());
         reductionManager->declareProperty(algProp);
       } else {
         g_log.error() << "ERROR: Sensitivity beam center determination was required"
             " but no file was provided" << std::endl;
       }
    }

    effAlg->setPropertyValue("OutputSensitivityWorkspace", outputSensitivityWS);
    effAlg->setPropertyValue("ReductionProperties", reductionManagerName);

    AlgorithmProperty *algProp = new AlgorithmProperty("SensitivityAlgorithm");
    algProp->setValue(effAlg->toString());
    reductionManager->declareProperty(algProp);
  }
}
void SetupEQSANSReduction::setupTransmission(boost::shared_ptr<PropertyManager> reductionManager)
{
  const std::string reductionManagerName = getProperty("ReductionProperties");
  // Transmission options
  const bool thetaDependentTrans = getProperty("ThetaDependentTransmission");
  const std::string transMethod = getProperty("TransmissionMethod");
  const std::string darkCurrent = getPropertyValue("TransmissionDarkCurrentFile");
  const bool useSampleDC = getProperty("TransmissionUseSampleDC");

  // Transmission is entered by hand
  if (boost::iequals(transMethod, "Value"))
  {
    const double transValue = getProperty("TransmissionValue");
    const double transError = getProperty("TransmissionError");
    if (!isEmpty(transValue) && !isEmpty(transError))
    {
      IAlgorithm_sptr transAlg = createChildAlgorithm("ApplyTransmissionCorrection");
      transAlg->setProperty("TransmissionValue", transValue);
      transAlg->setProperty("TransmissionError", transError);
      transAlg->setProperty("ThetaDependent", thetaDependentTrans);

      AlgorithmProperty *algProp = new AlgorithmProperty("TransmissionAlgorithm");
      algProp->setValue(transAlg->toString());
      reductionManager->declareProperty(algProp);
    } else {
      g_log.information("SetupEQSANSReduction [TransmissionAlgorithm]:"
          "expected transmission/error values and got empty values");
    }
  }
  // Direct beam method for transmission determination
  else if (boost::iequals(transMethod, "DirectBeam"))
  {
    const std::string sampleFilename = getPropertyValue("TransmissionSampleDataFile");
    const std::string emptyFilename = getPropertyValue("TransmissionEmptyDataFile");
    const double beamRadius = getProperty("TransmissionBeamRadius");
    const bool fitFramesTogether = getProperty("FitFramesTogether");
    const double beamX = getProperty("TransmissionBeamCenterX");
    const double beamY = getProperty("TransmissionBeamCenterY");
    const std::string centerMethod = getPropertyValue("TransmissionBeamCenterMethod");

    IAlgorithm_sptr transAlg = createChildAlgorithm("EQSANSDirectBeamTransmission");
    transAlg->setProperty("FitFramesTogether", fitFramesTogether);
    transAlg->setProperty("SampleDataFilename", sampleFilename);
    transAlg->setProperty("EmptyDataFilename", emptyFilename);
    transAlg->setProperty("BeamRadius", beamRadius);
    transAlg->setProperty("DarkCurrentFilename", darkCurrent);
    transAlg->setProperty("UseSampleDarkCurrent", useSampleDC);

    // Beam center option for transmission data
    if (boost::iequals(centerMethod, "Value") && !isEmpty(beamX) && !isEmpty(beamY))
    {
      transAlg->setProperty("BeamCenterX", beamX);
      transAlg->setProperty("BeamCenterY", beamY);
    }
    else if (boost::iequals(centerMethod, "DirectBeam"))
    {
      const std::string beamCenterFile = getProperty("TransmissionBeamCenterFile");
      if (beamCenterFile.size()>0)
       {
         IAlgorithm_sptr ctrAlg = createChildAlgorithm("SANSBeamFinder");
         ctrAlg->setProperty("Filename", beamCenterFile);
         ctrAlg->setProperty("UseDirectBeamMethod", true);
         ctrAlg->setProperty("PersistentCorrection", false);
         ctrAlg->setPropertyValue("ReductionProperties", reductionManagerName);

         AlgorithmProperty *algProp = new AlgorithmProperty("TransmissionBeamCenterAlgorithm");
         algProp->setValue(ctrAlg->toString());
         reductionManager->declareProperty(algProp);
       } else {
         g_log.error() << "ERROR: Transmission beam center determination was required"
             " but no file was provided" << std::endl;
       }
    }
    transAlg->setProperty("ThetaDependent", thetaDependentTrans);
    AlgorithmProperty *algProp = new AlgorithmProperty("TransmissionAlgorithm");
    algProp->setValue(transAlg->toString());
    reductionManager->declareProperty(algProp);
  }
}

void SetupEQSANSReduction::setupBackground(boost::shared_ptr<PropertyManager> reductionManager)
{
  const std::string reductionManagerName = getProperty("ReductionProperties");
  // Background
  const std::string backgroundFile = getPropertyValue("BackgroundFiles");
  if (backgroundFile.size() > 0)
    reductionManager->declareProperty(new PropertyWithValue<std::string>("BackgroundFiles", backgroundFile) );
  else
    return;

  const std::string darkCurrent = getPropertyValue("BckTransmissionDarkCurrentFile");
  const bool bckThetaDependentTrans = getProperty("BckThetaDependentTransmission");
  const std::string bckTransMethod = getProperty("BckTransmissionMethod");
  if (boost::iequals(bckTransMethod, "Value"))
  {
    const double transValue = getProperty("BckTransmissionValue");
    const double transError = getProperty("BckTransmissionError");
    if (!isEmpty(transValue) && !isEmpty(transError))
    {
      IAlgorithm_sptr transAlg = createChildAlgorithm("ApplyTransmissionCorrection");
      transAlg->setProperty("TransmissionValue", transValue);
      transAlg->setProperty("TransmissionError", transError);
      transAlg->setProperty("ThetaDependent", bckThetaDependentTrans);

      AlgorithmProperty *algProp = new AlgorithmProperty("BckTransmissionAlgorithm");
      algProp->setValue(transAlg->toString());
      reductionManager->declareProperty(algProp);
    } else {
      g_log.information("SetupEQSANSReduction [BckTransmissionAlgorithm]: "
          "expected transmission/error values and got empty values");
    }
  }
  else if (boost::iequals(bckTransMethod, "DirectBeam"))
  {
    const std::string sampleFilename = getPropertyValue("BckTransmissionSampleDataFile");
    const std::string emptyFilename = getPropertyValue("BckTransmissionEmptyDataFile");
    const double beamRadius = getProperty("BckTransmissionBeamRadius");
    const double beamX = getProperty("BckTransmissionBeamCenterX");
    const double beamY = getProperty("BckTransmissionBeamCenterY");
    const bool thetaDependentTrans = getProperty("BckThetaDependentTransmission");
    const bool useSampleDC = getProperty("TransmissionUseSampleDC");
    const bool fitFramesTogether = getProperty("BckFitFramesTogether");

    IAlgorithm_sptr transAlg = createChildAlgorithm("EQSANSDirectBeamTransmission");
    transAlg->setProperty("FitFramesTogether", fitFramesTogether);
    transAlg->setProperty("SampleDataFilename", sampleFilename);
    transAlg->setProperty("EmptyDataFilename", emptyFilename);
    transAlg->setProperty("BeamRadius", beamRadius);
    transAlg->setProperty("DarkCurrentFilename", darkCurrent);
    transAlg->setProperty("UseSampleDarkCurrent", useSampleDC);

    // Beam center option for transmission data
    const std::string centerMethod = getPropertyValue("BckTransmissionBeamCenterMethod");
    if (boost::iequals(centerMethod, "Value") && !isEmpty(beamX) && !isEmpty(beamY))
    {
      transAlg->setProperty("BeamCenterX", beamX);
      transAlg->setProperty("BeamCenterY", beamY);
    }
    else if (boost::iequals(centerMethod, "DirectBeam"))
    {
      const std::string beamCenterFile = getProperty("BckTransmissionBeamCenterFile");
      if (beamCenterFile.size()>0)
       {
         IAlgorithm_sptr ctrAlg = createChildAlgorithm("SANSBeamFinder");
         ctrAlg->setProperty("Filename", beamCenterFile);
         ctrAlg->setProperty("UseDirectBeamMethod", true);
         ctrAlg->setProperty("PersistentCorrection", false);
         ctrAlg->setPropertyValue("ReductionProperties", reductionManagerName);

         AlgorithmProperty *algProp = new AlgorithmProperty("BckTransmissionBeamCenterAlgorithm");
         algProp->setValue(ctrAlg->toString());
         reductionManager->declareProperty(algProp);
       } else {
         g_log.error() << "ERROR: Beam center determination was required"
             " but no file was provided" << std::endl;
       }
    }
    transAlg->setProperty("DarkCurrentFilename", darkCurrent);
    transAlg->setProperty("ThetaDependent", thetaDependentTrans);
    AlgorithmProperty *algProp = new AlgorithmProperty("BckTransmissionAlgorithm");
    algProp->setValue(transAlg->toString());
    reductionManager->declareProperty(algProp);
  }

}
} // namespace WorkflowAlgorithms
} // namespace Mantid

