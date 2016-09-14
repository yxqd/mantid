#ifndef ALGORITHMS_TEST_REFLECTOMETRYREDUCTIONONETEST_H_
#define ALGORITHMS_TEST_REFLECTOMETRYREDUCTIONONETEST_H_

#include <cxxtest/TestSuite.h>
#include <algorithm>
#include "MantidAlgorithms/ReflectometryReductionOne.h"
#include "MantidAPI/AlgorithmManager.h"
#include "MantidAPI/Axis.h"
#include "MantidAPI/FrameworkManager.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidTestHelpers/WorkspaceCreationHelper.h"
#include "MantidGeometry/Instrument/ReferenceFrame.h"
#include "MantidGeometry/Instrument.h"
#include <algorithm>

using namespace Mantid;
using namespace Mantid::Kernel;
using namespace Mantid::API;
using namespace Mantid::Algorithms;
using namespace Mantid::Geometry;
using namespace WorkspaceCreationHelper;

class ReflectometryReductionOneTest : public CxxTest::TestSuite {
private:
  MatrixWorkspace_sptr m_pointDetectorWS;
  MatrixWorkspace_sptr m_multiDetectorWS;
  MatrixWorkspace_sptr m_wavelengthWS;

public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static ReflectometryReductionOneTest *createSuite() {
    return new ReflectometryReductionOneTest();
  }
  static void destroySuite(ReflectometryReductionOneTest *suite) {
    delete suite;
  }

  ReflectometryReductionOneTest() {
    FrameworkManager::Instance();
    // A point detector ws 
    m_pointDetectorWS = create2DWorkspaceWithReflectometryInstrument();
    // A multi detector ws
    m_multiDetectorWS =
        create2DWorkspaceWithReflectometryInstrumentMultiDetector();
    // A workspace in wavelength
    m_wavelengthWS = Create2DWorkspace154(1, 10, true);
    m_wavelengthWS->setInstrument(m_pointDetectorWS->getInstrument());
    m_wavelengthWS->getAxis(0)->setUnit("Wavelength");

  }

  IAlgorithm_sptr construct_standard_algorithm() {
    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setRethrows(true);
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_pointDetectorWS);
    alg->setProperty("WavelengthMin", 1.0);
    alg->setProperty("WavelengthMax", 2.0);
    alg->setProperty("I0MonitorIndex", 1);
    alg->setProperty("MonitorBackgroundWavelengthMin", 1.0);
    alg->setProperty("MonitorBackgroundWavelengthMax", 2.0);
    alg->setProperty("MonitorIntegrationWavelengthMin", 1.2);
    alg->setProperty("MonitorIntegrationWavelengthMax", 1.5);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setPropertyValue("ProcessingInstructions", "1");
    alg->setPropertyValue("OutputWorkspace", "x");
    alg->setPropertyValue("OutputWorkspaceWavelength", "y");
    alg->setRethrows(true);
    return alg;
  }

  void test_execute() {
    auto alg = construct_standard_algorithm();
    TS_ASSERT_THROWS_NOTHING(alg->execute());
    MatrixWorkspace_sptr workspaceInQ = alg->getProperty("OutputWorkspace");
    MatrixWorkspace_sptr workspaceInLam =
        alg->getProperty("OutputWorkspaceWavelength");
    const double theta = alg->getProperty("ThetaOut");
    UNUSED_ARG(theta)
    UNUSED_ARG(workspaceInQ)
    UNUSED_ARG(workspaceInLam)
  }

  /// Conversion to wavelength

  void test_tolam() {
    MatrixWorkspace_sptr toConvert = m_pointDetectorWS;
    std::vector<int> detectorIndexRange;
    size_t workspaceIndexToKeep1 = 1;
    const int monitorIndex = 0;

    specnum_t specId1 =
        toConvert->getSpectrum(workspaceIndexToKeep1).getSpectrumNo();
    specnum_t monitorSpecId =
        toConvert->getSpectrum(monitorIndex).getSpectrumNo();

    // Define one spectra to keep
    detectorIndexRange.push_back(static_cast<int>(workspaceIndexToKeep1));
    std::stringstream buffer;
    buffer << workspaceIndexToKeep1;
    const std::string detectorIndexRangesStr = buffer.str();

    // Define a wavelength range for the detector workspace
    const double wavelengthMin = 1.0;
    const double wavelengthMax = 15;
    const double backgroundWavelengthMin = 17;
    const double backgroundWavelengthMax = 20;

    ReflectometryReductionOne alg;

    // Run the conversion.
    ReflectometryWorkflowBase::DetectorMonitorWorkspacePair inLam =
        alg.toLam(toConvert, detectorIndexRangesStr, monitorIndex,
                  boost::tuple<double, double>(wavelengthMin, wavelengthMax),
                  boost::tuple<double, double>(backgroundWavelengthMin,
                                               backgroundWavelengthMax));

    // Unpack the results
    MatrixWorkspace_sptr detectorWS = inLam.get<0>();
    MatrixWorkspace_sptr monitorWS = inLam.get<1>();

    /* ---------- Checks for the detector workspace ------------------*/

    // Check units.
    TS_ASSERT_EQUALS("Wavelength", detectorWS->getAxis(0)->unit()->unitID());

    // Check the number of spectrum kept.
    TS_ASSERT_EQUALS(1, detectorWS->getNumberHistograms());

    auto map = detectorWS->getSpectrumToWorkspaceIndexMap();
    // Check the spectrum Nos retained.
    TS_ASSERT_EQUALS(map[specId1], 0);

    // Check the cropped x range
    Mantid::MantidVec copyX = detectorWS->readX(0);
    std::sort(copyX.begin(), copyX.end());
    TS_ASSERT(copyX.front() >= wavelengthMin);
    TS_ASSERT(copyX.back() <= wavelengthMax);

    /* ------------- Checks for the monitor workspace --------------------*/
    // Check units.
    TS_ASSERT_EQUALS("Wavelength", monitorWS->getAxis(0)->unit()->unitID());

    // Check the number of spectrum kept. This should only ever be 1.
    TS_ASSERT_EQUALS(1, monitorWS->getNumberHistograms());

    map = monitorWS->getSpectrumToWorkspaceIndexMap();
    // Check the spectrum Nos retained.
    TS_ASSERT_EQUALS(map[monitorSpecId], 0);
  }

  void test_wavelength_conversion_1() {
    // Monitors : No
    // Direct beam normalization : No
    // Processing instructions : 1
    // Analysis mode : multi detector

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_multiDetectorWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "MultiDetectorAnalysis");
    alg->setProperty("NormalizeByIntegratedMonitors", "0");
    alg->setPropertyValue("ProcessingInstructions", "1");
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");
    alg->execute();
    // We are only interested in workspace in wavelength
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");

    TS_ASSERT(outLam);
    TS_ASSERT_EQUALS(outLam->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outLam->blocksize(), 8);
    TS_ASSERT(outLam->readX(0)[0] >= 1.5);
    TS_ASSERT(outLam->readX(0)[7] <= 15.0);
    TS_ASSERT_DELTA(outLam->readY(0)[0], 3.1530, 0.0001);
    TS_ASSERT_DELTA(outLam->readY(0)[7], 3.1530, 0.0001);
  }

  void test_wavelength_conversion_2() {
    // Monitors : No
    // Direct beam normalization : No
    // Processing instructions : 1+2
    // Analysis mode : multi detector

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_multiDetectorWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "MultiDetectorAnalysis");
    alg->setProperty("NormalizeByIntegratedMonitors", "0");
    alg->setPropertyValue("ProcessingInstructions", "1+2");
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");
    alg->execute();
    // We are only interested in workspace in wavelength
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");

    TS_ASSERT_EQUALS(outLam->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outLam->blocksize(), 8);
    TS_ASSERT(outLam->readX(0)[0] >= 1.5);
    TS_ASSERT(outLam->readX(0)[7] <= 15.0);
    // Y counts, should be 3.1530 * 2 (see test_wavelength_conversion_1())
    TS_ASSERT_DELTA(outLam->readY(0)[0], 6.3060, 0.0001);
    TS_ASSERT_DELTA(outLam->readY(0)[7], 6.3060, 0.0001);
  }

  void test_wavelength_conversion_3() {
    // Monitors : No
    // Direct beam normalization : No
    // Processing instructions : 1-3 (i.e. all detectors)
    // Analysis mode : multi detector

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_multiDetectorWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "MultiDetectorAnalysis");
    alg->setProperty("NormalizeByIntegratedMonitors", "0");
    alg->setPropertyValue("ProcessingInstructions", "1-3");
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");
    alg->execute();
    // We are only interested in workspace in wavelength
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");

    TS_ASSERT_EQUALS(outLam->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outLam->blocksize(), 8);
    TS_ASSERT(outLam->readX(0)[0] >= 1.5);
    TS_ASSERT(outLam->readX(0)[7] <= 15.0);
    // Y counts, should be 3.1530 * 3 (see test_wavelength_conversion_1())
    TS_ASSERT_DELTA(outLam->readY(0)[0], 9.4590, 0.0001);
    TS_ASSERT_DELTA(outLam->readY(0)[7], 9.4590, 0.0001);
  }

  void test_wavelength_conversion_4() {
    // Monitors : No
    // Direct beam normalization : Yes
    // Processing instructions : 1
    // Region of direct beam : 2-3
    // Analysis mode : multi detector

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_multiDetectorWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "MultiDetectorAnalysis");
    alg->setProperty("NormalizeByIntegratedMonitors", "0");
    alg->setPropertyValue("ProcessingInstructions", "1");
    alg->setPropertyValue("RegionOfDirectBeam", "2-3");
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");
    alg->execute();
    // We are only interested in workspace in wavelength
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");

    TS_ASSERT_EQUALS(outLam->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outLam->blocksize(), 8);
    TS_ASSERT(outLam->readX(0)[0] >= 1.5);
    TS_ASSERT(outLam->readX(0)[7] <= 15.0);
    // Y counts, should be 0.5 = 1 (from detector ws) / 2 (from direct beam)
    TS_ASSERT_DELTA(outLam->readY(0)[0], 0.5, 0.0001);
    TS_ASSERT_DELTA(outLam->readY(0)[7], 0.5, 0.0001);
  }

  void test_wavelength_conversion_5() {
    // Monitors : Yes (0)
    // MonitorBackgroundWavelengthMin : Not given
    // MonitorBackgroundWavelengthMax : Not given
    // Normalize by integrated monitors : No
    // Direct beam normalization : No
    // Processing instructions : 1
    // Region of direct beam : No
    // Analysis mode : multi detector

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_multiDetectorWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "MultiDetectorAnalysis");
    alg->setProperty("I0MonitorIndex", "0");
    alg->setProperty("NormalizeByIntegratedMonitors", "0");
    alg->setPropertyValue("ProcessingInstructions", "1");
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");
    alg->execute();
    // We are only interested in workspace in wavelength
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");

    TS_ASSERT_EQUALS(outLam->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outLam->blocksize(), 8);
    TS_ASSERT(outLam->readX(0)[0] >= 1.5);
    TS_ASSERT(outLam->readX(0)[7] <= 15.0);
    // No monitors considered because MonitorBackgroundWavelengthMin
    // and MonitorBackgroundWavelengthMax were not set
    // even if I0MonitorIndex was given
    // Y counts must be 3.1530
    TS_ASSERT_DELTA(outLam->readY(0)[0], 3.1530, 0.0001);
    TS_ASSERT_DELTA(outLam->readY(0)[7], 3.1530, 0.0001);
  }

  void test_wavelength_conversion_6() {
    // Monitors : Yes (0)
    // MonitorBackgroundWavelengthMin : Yes
    // MonitorBackgroundWavelengthMax : Not given
    // Normalize by integrated monitors : No
    // Direct beam normalization : No
    // Processing instructions : 1
    // Region of direct beam : No
    // Analysis mode : multi detector

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_multiDetectorWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "MultiDetectorAnalysis");
    alg->setProperty("I0MonitorIndex", "0");
    alg->setProperty("MonitorBackgroundWavelengthMin", 1.5);
    alg->setProperty("NormalizeByIntegratedMonitors", "0");
    alg->setPropertyValue("ProcessingInstructions", "1");
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");
    alg->execute();
    // We are only interested in workspace in wavelength
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");

    TS_ASSERT_EQUALS(outLam->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outLam->blocksize(), 8);
    TS_ASSERT(outLam->readX(0)[0] >= 1.5);
    TS_ASSERT(outLam->readX(0)[7] <= 15.0);
    // No monitors considered because MonitorBackgroundWavelengthMax
    // was not set
    // even if I0MonitorIndex was given
    // Y counts must be 3.1530
    TS_ASSERT_DELTA(outLam->readY(0)[0], 3.1530, 0.0001);
    TS_ASSERT_DELTA(outLam->readY(0)[7], 3.1530, 0.0001);
  }

  void test_wavelength_conversion_7() {
    // Monitors : Yes (0)
    // MonitorBackgroundWavelengthMin : Not given
    // MonitorBackgroundWavelengthMax : 15.0
    // Normalize by integrated monitors : No
    // Direct beam normalization : No
    // Processing instructions : 1
    // Region of direct beam : No
    // Analysis mode : multi detector

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_multiDetectorWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "MultiDetectorAnalysis");
    alg->setProperty("I0MonitorIndex", "0");
    alg->setProperty("MonitorBackgroundWavelengthMax", 15.0);
    alg->setProperty("NormalizeByIntegratedMonitors", "0");
    alg->setPropertyValue("ProcessingInstructions", "1");
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");
    alg->execute();
    // We are only interested in workspace in wavelength
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");

    TS_ASSERT_EQUALS(outLam->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outLam->blocksize(), 8);
    TS_ASSERT(outLam->readX(0)[0] >= 1.5);
    TS_ASSERT(outLam->readX(0)[7] <= 15.0);
    // No monitors considered because MonitorBackgroundWavelengthMin
    // was not set
    // even if I0MonitorIndex was given
    // Y counts must be 3.1530
    TS_ASSERT_DELTA(outLam->readY(0)[0], 3.1530, 0.0001);
    TS_ASSERT_DELTA(outLam->readY(0)[7], 3.1530, 0.0001);
  }

  void test_wavelength_conversion_8() {
    // Monitors : Yes (0)
    // MonitorBackgroundWavelengthMin : 0.5
    // MonitorBackgroundWavelengthMax : 3.0
    // Normalize by integrated monitors : No
    // Direct beam normalization : No
    // Processing instructions : 1
    // Region of direct beam : No
    // Analysis mode : multi detector

    // Modify counts in monitor (only for this test)
    // Modify counts only for range that will be fitted
    auto inputWS = m_multiDetectorWS;
    MantidVec &Y = inputWS->dataY(0);
    std::fill(Y.begin(), Y.begin()+2, 1.0);

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", inputWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "MultiDetectorAnalysis");
    alg->setProperty("I0MonitorIndex", "0");
    alg->setProperty("MonitorBackgroundWavelengthMin", 0.5);
    alg->setProperty("MonitorBackgroundWavelengthMax", 3.0);
    alg->setProperty("NormalizeByIntegratedMonitors", "0");
    alg->setPropertyValue("ProcessingInstructions", "1");
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");
    alg->execute();
    // We are only interested in workspace in wavelength
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");

    TS_ASSERT_EQUALS(outLam->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outLam->blocksize(), 8);
    TS_ASSERT(outLam->readX(0)[0] >= 1.5);
    TS_ASSERT(outLam->readX(0)[7] <= 15.0);
    // Expected values are 2.4996 = 3.15301 (detectors) / 1.26139 (monitors)
    TS_ASSERT_DELTA(outLam->readY(0)[0], 2.4996, 0.0001);
    TS_ASSERT_DELTA(outLam->readY(0)[7], 2.4996, 0.0001);
  }

  void test_wavelength_conversion_9() {
    // Monitors : Yes (0)
    // MonitorBackgroundWavelengthMin : 0.5
    // MonitorBackgroundWavelengthMax : 3.0
    // Normalize by integrated monitors : Yes
    // Direct beam normalization : No
    // Processing instructions : 1
    // Region of direct beam : No
    // Analysis mode : multi detector

    // Modify counts in monitor (only for this test)
    // Modify counts only for range that will be fitted
    auto inputWS = m_multiDetectorWS;
    MantidVec &Y = inputWS->dataY(0);
    std::fill(Y.begin(), Y.begin()+2, 1.0);

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", inputWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "MultiDetectorAnalysis");
    alg->setProperty("I0MonitorIndex", "0");
    alg->setProperty("MonitorBackgroundWavelengthMin", 0.5);
    alg->setProperty("MonitorBackgroundWavelengthMax", 3.0);
    alg->setProperty("NormalizeByIntegratedMonitors", "1");
    alg->setPropertyValue("ProcessingInstructions", "1");
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");
    alg->execute();
    // We are only interested in workspace in wavelength
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");

    TS_ASSERT_EQUALS(outLam->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outLam->blocksize(), 8);
    TS_ASSERT(outLam->readX(0)[0] >= 1.5);
    TS_ASSERT(outLam->readX(0)[7] <= 15.0);
    // Expected values are 0.3124 = 3.15301 (detectors) / (1.26139*8) (monitors)
    TS_ASSERT_DELTA(outLam->readY(0)[0], 0.3124, 0.0001);
    TS_ASSERT_DELTA(outLam->readY(0)[7], 0.3124, 0.0001);
  }

  void test_wavelength_conversion_10() {
    // Analysis mode : point
    // Monitors : No
    // Direct beam normalization : No
    // Processing instructions : 1+2 and 1

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_pointDetectorWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "PointDetectorAnalysis");
    alg->setProperty("NormalizeByIntegratedMonitors", "0");
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");

    // Must throw as spectrum 2 is not defined
    alg->setPropertyValue("ProcessingInstructions", "1+2");
    TS_ASSERT_THROWS_ANYTHING(alg->execute());
    // Fine, we are summing first spectrum twice
    alg->setPropertyValue("ProcessingInstructions", "1+1");
    TS_ASSERT_THROWS_NOTHING(alg->execute());
  }

  void test_wavelength_conversion_11() {
    // Analysis mode : point
    // Monitors : No
    // Direct beam normalization : 1-1, 1-2
    // Processing instructions : 1

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_pointDetectorWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "PointDetectorAnalysis");
    alg->setProperty("NormalizeByIntegratedMonitors", "0");
    alg->setPropertyValue("ProcessingInstructions", "1");
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");

    // Must throw region of interest is incompatible with point detector
    alg->setPropertyValue("RegionOfDirectBeam", "1-2");
    TS_ASSERT_THROWS_ANYTHING(alg->execute());
    // Even if we set only the existing spectrum
    alg->setPropertyValue("RegionOfDirectBeam", "1-1");
    TS_ASSERT_THROWS_ANYTHING(alg->execute());
  }

  /// Transmission correction

  void test_transmission_correction_run() {
    // CorrectDetectorPositions: False

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_wavelengthWS);
    alg->setProperty("FirstTransmissionRun", m_wavelengthWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "PointDetectorAnalysis");
    alg->setProperty("ProcessingInstructions", "100");
    alg->setProperty("CorrectDetectorPositions", false);
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");
    alg->execute();
    // We are only interested in workspace in wavelength
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");

    TS_ASSERT_EQUALS(outLam->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outLam->blocksize(), 10);
    // Expected values are 1 = m_wavelength / m_wavelength
    TS_ASSERT_DELTA(outLam->readY(0)[0], 1.0000, 0.0001);
    TS_ASSERT_DELTA(outLam->readY(0)[7], 1.0000, 0.0001);
  }

  void test_transmission_correction_exponential() {
    // CorrectionAlgorithm: ExponentialCorrection
    
    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_wavelengthWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "PointDetectorAnalysis");
    alg->setProperty("ProcessingInstructions", "100");
    alg->setProperty("CorrectDetectorPositions", false);
    alg->setProperty("CorrectionAlgorithm", "ExponentialCorrection");
    alg->setProperty("C0", 0.2);
    alg->setProperty("C1", 0.1);
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");
    alg->execute();
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");

    TS_ASSERT_EQUALS(outLam->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outLam->blocksize(), 10);
    // Expected values are 29.0459 and 58.4912
    TS_ASSERT_DELTA(outLam->readY(0)[0], 29.0459, 0.0001);
    TS_ASSERT_DELTA(outLam->readY(0)[7], 58.4912, 0.0001);
  }

  void test_transmission_correction_polynomial() {
    // CorrectionAlgorithm: PolynomialCorrection

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_wavelengthWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "PointDetectorAnalysis");
    alg->setProperty("ProcessingInstructions", "100");
    alg->setProperty("CorrectDetectorPositions", false);
    alg->setProperty("CorrectionAlgorithm", "PolynomialCorrection");
    alg->setProperty("Polynomial", "0.1,0.3,0.5");
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");
    alg->execute();
    // We are only interested in workspace in wavelength
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");

    TS_ASSERT_EQUALS(outLam->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outLam->blocksize(), 10);
    // Expected values are 2.9851 and 0.1289
    TS_ASSERT_DELTA(outLam->readY(0)[0], 2.9851, 0.0001);
    TS_ASSERT_DELTA(outLam->readY(0)[7], 0.1289, 0.0001);
  }

  /// Detector position correction

  void test_detector_position_correction_point_detector() {
    // Analysis mode: point detector

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_pointDetectorWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "PointDetectorAnalysis");
    alg->setProperty("NormalizeByIntegratedMonitors", "0");
    alg->setProperty("ProcessingInstructions", "1");
    alg->setProperty("CorrectDetectorPositions", true);
    alg->setProperty("ThetaIn", 1.5);
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");
    alg->execute();
    // We are only interested in workspace in wavelength
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");

    TS_ASSERT_EQUALS(outLam->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outLam->blocksize(), 23);
    // Check detector positions have changed from original and are correct
    TS_ASSERT_DIFFERS(outLam->getDetector(0)->getPos(), m_pointDetectorWS->getDetector(0)->getPos());
    TS_ASSERT_EQUALS(outLam->getDetector(0)->getPos(), V3D(14, 0, 0));
    // Expected values are 2.0000 = 2.0000 (detectors) / 1.0000 (monitors)
    TS_ASSERT_DELTA(outLam->readY(0)[0], 2.0000, 0.0001);
    TS_ASSERT_DELTA(outLam->readY(0)[7], 2.0000, 0.0001);
  }

  void test_detector_position_correction_multi_detector() {
    // Analysis mode: multi detector

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_multiDetectorWS);
    alg->setProperty("WavelengthMin", 1.5);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("AnalysisMode", "MultiDetectorAnalysis");
    alg->setProperty("NormalizeByIntegratedMonitors", "0");
    alg->setProperty("ProcessingInstructions", "1");
    alg->setProperty("CorrectDetectorPositions", true);
    alg->setProperty("ThetaIn", 1.5);
    alg->setPropertyValue("OutputWorkspace", "IvsQ");
    alg->setPropertyValue("OutputWorkspaceWavelength", "IvsLam");
    alg->execute();
    // We are only interested in workspace in wavelength
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");

    TS_ASSERT_EQUALS(outLam->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outLam->blocksize(), 8);
    // Check detector positions have changed from original and are correct
    TS_ASSERT_DIFFERS(outLam->getDetector(0)->getPos(),
                      m_multiDetectorWS->getDetector(0)->getPos());
    TS_ASSERT_EQUALS(outLam->getDetector(0)->getPos(), V3D(20, 0.13093, 0));
    // Expected values are 3.1530 = 3.15301 (detectors) / 1.0000 (monitors)
    TS_ASSERT_DELTA(outLam->readY(0)[0], 3.1530, 0.0001);
    TS_ASSERT_DELTA(outLam->readY(0)[7], 3.1530, 0.0001);
  }

  /// Calculate theta

  void test_calculate_theta() {

    auto alg = construct_standard_algorithm();

    alg->execute();
    // Should not throw

    const double outTwoTheta = alg->getProperty("ThetaOut");
    TS_ASSERT_DELTA(45.0, outTwoTheta, 0.00001);
  }

  /// Conversion to momentum transfer

  void test_source_rotation_after_second_reduction() {
    // set up the axis for the instrument
    Instrument_sptr instrument = boost::make_shared<Instrument>();
    instrument->setReferenceFrame(boost::make_shared<ReferenceFrame>(
        Y /*up*/, Z /*along*/, Right, "0,0,0"));

    // add a source
    ObjComponent *source = new ObjComponent("source");
    source->setPos(V3D(0, 0, -1));
    instrument->add(source);
    instrument->markAsSource(source);

    // add a sample
    ObjComponent *sample = new ObjComponent("some-surface-holder");
    sample->setPos(V3D(0, 0, 0));
    instrument->add(sample);
    instrument->markAsSamplePos(sample);

    // add a detector
    Detector *det = new Detector("point-detector", 1, NULL);
    det->setPos(V3D(0, 1, 1));
    instrument->add(det);
    instrument->markAsDetector(det);

    // set the instrument to this workspace
    m_pointDetectorWS->setInstrument(instrument);
    // set this detector ready for processing instructions
    m_pointDetectorWS->getSpectrum(0).setDetectorID(det->getID());

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setRethrows(true);
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_pointDetectorWS);
    alg->setProperty("WavelengthMin", 1.0);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("I0MonitorIndex", 0);
    alg->setProperty("MonitorBackgroundWavelengthMin", 0.0);
    alg->setProperty("MonitorBackgroundWavelengthMax", 0.0);
    alg->setProperty("MonitorIntegrationWavelengthMin", 0.0);
    alg->setProperty("MonitorIntegrationWavelengthMax", 0.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("NormalizeByIntegratedMonitors", false);
    alg->setProperty("CorrectDetectorPositions", true);
    alg->setProperty("CorrectionAlgorithm", "None");
    alg->setPropertyValue("ProcessingInstructions", "1");
    alg->setPropertyValue("OutputWorkspace", "x");
    alg->setPropertyValue("OutputWorkspaceWavelength", "y");
    alg->setRethrows(true);
    TS_ASSERT_THROWS_NOTHING(alg->execute());
    MatrixWorkspace_sptr outLam = alg->getProperty("OutputWorkspaceWavelength");
    MatrixWorkspace_sptr outQ = alg->getProperty("OutputWorkspace");

    TS_ASSERT_EQUALS(m_pointDetectorWS->getInstrument()->getSource()->getPos(),
                     outLam->getInstrument()->getSource()->getPos());
    TS_ASSERT_EQUALS(outLam->getInstrument()->getSource()->getPos(),
                     outQ->getInstrument()->getSource()->getPos());
  }

  void test_scale_step() {
    auto alg = construct_standard_algorithm();
    auto inWS =
        WorkspaceCreationHelper::create2DWorkspaceWithReflectometryInstrument(
            2.0);
    inWS->getAxis(0)->setUnit("Wavelength");
    alg->setProperty("InputWorkspace", inWS);
    alg->setProperty("ScaleFactor", 1.0);
    alg->setProperty("ThetaIn", 1.5);
    alg->setProperty("OutputWorkspace", "Test");
    TS_ASSERT_THROWS_NOTHING(alg->execute());
    MatrixWorkspace_sptr nonScaledWS = alg->getProperty("OutputWorkspace");
    alg->setProperty("InputWorkspace", inWS);
    alg->setProperty("ScaleFactor", 0.5);
    alg->setProperty("OutputWorkspace", "scaledTest");
    TS_ASSERT_THROWS_NOTHING(alg->execute());
    MatrixWorkspace_sptr scaledWS = alg->getProperty("OutputWorkspace");
    // compare y data instead of workspaces.
    auto scaledYData = scaledWS->readY(0);
    auto nonScaledYData = nonScaledWS->readY(0);
    TS_ASSERT_EQUALS(scaledYData.front(), 2 * nonScaledYData.front());
    TS_ASSERT_EQUALS(scaledYData[scaledYData.size() / 2],
                     2 * nonScaledYData[nonScaledYData.size() / 2]);
    TS_ASSERT_EQUALS(scaledYData.back(), 2 * nonScaledYData.back());
    // Remove workspace from the data service.
    AnalysisDataService::Instance().remove("Test");
    AnalysisDataService::Instance().remove("scaledTest");
  }

  void test_rebin_in_Q_params_not_provided() {
    auto alg = construct_standard_algorithm();
    auto inWS = m_wavelengthWS;
    // this instrument does not have a "slit-gap" property
    // defined in the IPF, so CalculateResolution should throw.
    // Setup bad bin edges, Rebin will throw (not CalculateResolution?)
    inWS->dataX(0).assign(inWS->readX(0).size(), inWS->readX(0)[0]);
    alg->setProperty("InputWorkspace", inWS);
    alg->setProperty("OutputWorkspace", "rebinnedWS");
    TS_ASSERT_THROWS(alg->execute(), std::invalid_argument);
  }

  void test_rebin_in_Q_partial_params_provided() {
    auto alg = construct_standard_algorithm();
    alg->setProperty("InputWorkspace", m_wavelengthWS);
    alg->setProperty("MomentumTransferMaximum", 15.0);
    alg->setProperty("OutputWorkspace", "rebinnedWS");
    TS_ASSERT_THROWS_NOTHING(alg->execute());
    MatrixWorkspace_sptr rebinnedIvsQWS = alg->getProperty("OutputWorkspace");
    auto xData = rebinnedIvsQWS->readX(0);
    // based off the equation for logarithmic binning X(i+1)=X(i)(1+|dX|)
    double binWidthFromLogarithmicEquation = fabs((xData[1] / xData[0]) - 1);
    TSM_ASSERT_DELTA("DQQ should be the same as abs(x[1]/x[0] - 1)",
                     binWidthFromLogarithmicEquation, 0.1, 1e-06);
    TSM_ASSERT_DELTA("Qmax should be the same as last Params entry (5.0)",
                     xData.back(), 15.0, 1e-06);
  }

  void test_rebin_in_Q_logarithmic_rebinning() {
    auto alg = construct_standard_algorithm();
    alg->setProperty("InputWorkspace", m_wavelengthWS);
    alg->setProperty("MomentumTransferMinimum", 1.0);
    alg->setProperty("MomentumTransferStep", 0.2);
    alg->setProperty("MomentumTransferMaximum", 5.0);
    alg->setProperty("OutputWorkspace", "rebinnedWS");
    TS_ASSERT_THROWS_NOTHING(alg->execute());
    MatrixWorkspace_sptr rebinnedIvsQWS = alg->getProperty("OutputWorkspace");
    auto xData = rebinnedIvsQWS->readX(0);
    TSM_ASSERT_EQUALS("QMin should be the same as first Param entry (1.0)",
                      xData[0], 1.0);
    // based off the equation for logarithmic binning X(i+1)=X(i)(1+|dX|)
    double binWidthFromLogarithmicEquation = fabs((xData[1] / xData[0]) - 1);
    TSM_ASSERT_DELTA("DQQ should be the same as abs(x[1]/x[0] - 1)",
                     binWidthFromLogarithmicEquation, 0.2, 1e-06);
    TSM_ASSERT_EQUALS("QMax should be the same as last Param entry",
                      xData.back(), 5.0);
  }

  void test_rebin_in_Q_linear_rebinning() {
    auto alg = construct_standard_algorithm();
    alg->setProperty("InputWorkspace", m_wavelengthWS);
    alg->setProperty("MomentumTransferMinimum", 1.577);
    alg->setProperty("MomentumTransferStep", -0.2);
    alg->setProperty("MomentumTransferMaximum", 5.233);
    alg->setProperty("OutputWorkspace", "rebinnedWS");
    TS_ASSERT_THROWS_NOTHING(alg->execute());
    MatrixWorkspace_sptr rebinnedIvsQWS = alg->getProperty("OutputWorkspace");
    auto xData = rebinnedIvsQWS->readX(0);
    TSM_ASSERT_DELTA("QMin should be the same as the first Param entry (1.577)",
                     xData[0], 1.577, 1e-06);
    TSM_ASSERT_DELTA("DQQ should the same as 0.2", xData[1] - xData[0], 0.2,
                     1e-06);
    TSM_ASSERT_DELTA("QMax should be the same as the last Param entry (5.233)",
                     xData.back(), 5.233, 1e-06);
  }

  void test_Q_range() {
    // set up the axis for the instrument
    Instrument_sptr instrument = boost::make_shared<Instrument>();
    instrument->setReferenceFrame(boost::make_shared<ReferenceFrame>(
        Y /*up*/, Z /*along*/, Right, "0,0,0"));

    // add a source
    ObjComponent *source = new ObjComponent("source");
    source->setPos(V3D(0, 0, -1));
    instrument->add(source);
    instrument->markAsSource(source);

    // add a sample
    ObjComponent *sample = new ObjComponent("some-surface-holder");
    sample->setPos(V3D(0, 0, 0));
    instrument->add(sample);
    instrument->markAsSamplePos(sample);

    // add a detector
    Detector *det = new Detector("point-detector", 1, NULL);
    det->setPos(V3D(0, 1, 1));
    instrument->add(det);
    instrument->markAsDetector(det);

    // set the instrument to this workspace
    m_pointDetectorWS->setInstrument(instrument);
    // set this detector ready for processing instructions
    m_pointDetectorWS->getSpectrum(0).setDetectorID(det->getID());

    auto alg = AlgorithmManager::Instance().create("ReflectometryReductionOne");
    alg->setRethrows(true);
    alg->setChild(true);
    alg->initialize();
    alg->setProperty("InputWorkspace", m_pointDetectorWS);
    alg->setProperty("WavelengthMin", 1.0);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("I0MonitorIndex", 0);
    alg->setProperty("MonitorBackgroundWavelengthMin", 0.0);
    alg->setProperty("MonitorBackgroundWavelengthMax", 0.0);
    alg->setProperty("MonitorIntegrationWavelengthMin", 0.0);
    alg->setProperty("MonitorIntegrationWavelengthMax", 0.0);
    alg->setProperty("MomentumTransferStep", 0.1);
    alg->setProperty("NormalizeByIntegratedMonitors", false);
    alg->setProperty("CorrectDetectorPositions", true);
    alg->setProperty("CorrectionAlgorithm", "None");
    alg->setPropertyValue("ProcessingInstructions", "1");
    alg->setPropertyValue("OutputWorkspace", "x");
    alg->setPropertyValue("OutputWorkspaceWavelength", "y");
    alg->setRethrows(true);
    TS_ASSERT_THROWS_NOTHING(alg->execute());

    // retrieve the IvsLam workspace
    MatrixWorkspace_sptr inLam = alg->getProperty("OutputWorkspaceWavelength");
    // retrieve the IvsQ workspace
    MatrixWorkspace_sptr inQ = alg->getProperty("OutputWorkspace");
    // retrieve our Theta
    double outTheta = alg->getProperty("ThetaOut");

    TS_ASSERT_DELTA(45.0, outTheta, 0.00001);
    TS_ASSERT_EQUALS(source->getPos(),
                     inQ->getInstrument()->getSource()->getPos());
    // convert from degrees to radians for sin() function
    double outThetaInRadians = outTheta * M_PI / 180;

    double lamMin = inLam->readX(0).front();
    double lamMax = inLam->readX(0).back();

    // Derive our QMin and QMax from the equation
    double qMinFromEQ = (4 * M_PI * sin(outThetaInRadians)) / lamMax;
    double qMaxFromEQ = (4 * M_PI * sin(outThetaInRadians)) / lamMin;

    // Get our QMin and QMax from the workspace
    auto qMinFromWS = inQ->readX(0).front();
    auto qMaxFromWS = inQ->readX(0).back();

    // Compare the two values (they should be identical)
    TS_ASSERT_DELTA(qMinFromEQ, qMinFromWS, 0.00001);
    TS_ASSERT_DELTA(qMaxFromEQ, qMaxFromWS, 0.00001);
  }
};

#endif /* ALGORITHMS_TEST_REFLECTOMETRYREDUCTIONONETEST_H_ */
