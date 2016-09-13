/*
 * CreateTransmissionWorkspaceTest.h
 *
 *  Created on: Jul 29, 2014
 *      Author: spu92482
 */

#ifndef ALGORITHMS_TEST_CREATETRANSMISSIONWORKSPACETEST_H_
#define ALGORITHMS_TEST_CREATETRANSMISSIONWORKSPACETEST_H_

#include <cxxtest/TestSuite.h>
#include <algorithm>
#include "MantidAlgorithms/ReflectometryReductionOne.h"
#include "MantidAPI/AlgorithmManager.h"
#include "MantidAPI/Axis.h"
#include "MantidAPI/FrameworkManager.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidTestHelpers/WorkspaceCreationHelper.h"
#include "MantidGeometry/Instrument/ReferenceFrame.h"

using namespace Mantid;
using namespace Mantid::Kernel;
using namespace Mantid::API;
using namespace Mantid::Algorithms;
using namespace WorkspaceCreationHelper;

class CreateTransmissionWorkspaceTest : public CxxTest::TestSuite {
private:
  MatrixWorkspace_sptr m_pointDetectorWS;
  MatrixWorkspace_sptr m_TOF;
  MatrixWorkspace_sptr m_NotTOF;

private:
  IAlgorithm_sptr construct_standard_algorithm() {
    auto alg =
        AlgorithmManager::Instance().create("CreateTransmissionWorkspaceAuto");
    alg->initialize();
    alg->setChild(true);
    alg->setProperty("FirstTransmissionRun", m_pointDetectorWS);
    alg->setProperty("WavelengthMin", 0.0);
    alg->setProperty("WavelengthMax", 1.0);
    alg->setProperty("I0MonitorIndex", 0);
    alg->setPropertyValue("ProcessingInstructions", "0, 1");
    alg->setProperty("MonitorBackgroundWavelengthMin", 0.0);
    alg->setProperty("MonitorBackgroundWavelengthMax", 1.0);
    alg->setProperty("MonitorIntegrationWavelengthMin", 0.0);
    alg->setProperty("MonitorIntegrationWavelengthMax", 1.0);
    alg->setPropertyValue("OutputWorkspace", "demo_ws");
    alg->setRethrows(true);
    return alg;
  }

public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static CreateTransmissionWorkspaceTest *createSuite() {
    return new CreateTransmissionWorkspaceTest();
  }
  static void destroySuite(CreateTransmissionWorkspaceTest *suite) {
    delete suite;
  }

  CreateTransmissionWorkspaceTest() {
    m_pointDetectorWS = create2DWorkspaceWithReflectometryInstrument();
    m_NotTOF = m_pointDetectorWS;
    m_NotTOF->getAxis(0)->setUnit("MomentumTransfer");
  }

  void test_check_first_transmission_workspace_not_tof_or_wavelength_throws() {
    auto alg = construct_standard_algorithm();
    TS_ASSERT_THROWS(alg->setProperty("FirstTransmissionRun", m_NotTOF),
                     std::invalid_argument);
  }

  void test_check_second_transmission_workspace_not_tof_throws() {
    auto alg = construct_standard_algorithm();
    TS_ASSERT_THROWS(alg->setProperty("SecondTransmissionRun", m_NotTOF),
                     std::invalid_argument);
  }

  void test_end_overlap_must_be_greater_than_start_overlap_or_throw() {
    auto alg = construct_standard_algorithm();
    alg->setProperty("FirstTransmissionRun", m_pointDetectorWS);
    alg->setProperty("SecondTransmissionRun", m_pointDetectorWS);
    MantidVec params = {0.0, 0.1, 1.0};
    alg->setProperty("Params", params);
    alg->setProperty("StartOverlap", 0.6);
    alg->setProperty("EndOverlap", 0.4);
    TS_ASSERT_THROWS(alg->execute(), std::invalid_argument);
  }

  void test_must_provide_wavelengths() {
    auto alg =
        AlgorithmManager::Instance().create("CreateTransmissionWorkspace");
    alg->initialize();
    alg->setChild(true);
    alg->setProperty("FirstTransmissionRun", m_pointDetectorWS);
    alg->setProperty("SecondTransmissionRun", m_pointDetectorWS);
    alg->setPropertyValue("OutputWorkspace", "demo_ws");
    alg->setRethrows(true);
    TS_ASSERT_THROWS(alg->execute(), std::runtime_error);

    alg->setProperty("FirstTransmissionRun", m_pointDetectorWS);
    alg->setProperty("SecondTransmissionRun", m_pointDetectorWS);
    alg->setProperty("WavelengthMin", 1.0);
    alg->setRethrows(true);
    TS_ASSERT_THROWS(alg->execute(), std::runtime_error);
  }

  void test_wavelength_min_greater_wavelength_max_throws() {
    auto alg = construct_standard_algorithm();
    alg->setProperty("WavelengthMin", 1.0);
    alg->setProperty("WavelengthMax", 0.0);
    TS_ASSERT_THROWS(alg->execute(), std::invalid_argument);
  }

  void
  test_monitor_background_wavelength_min_greater_monitor_background_wavelength_max_throws() {
    auto alg = construct_standard_algorithm();
    alg->setProperty("MonitorBackgroundWavelengthMin", 1.0);
    alg->setProperty("MonitorBackgroundWavelengthMax", 0.0);
    TS_ASSERT_THROWS(alg->execute(), std::invalid_argument);
  }

  void
  test_monitor_integration_wavelength_min_greater_monitor_integration_wavelength_max_throws() {
    auto alg = construct_standard_algorithm();
    alg->setProperty("MonitorIntegrationWavelengthMin", 1.0);
    alg->setProperty("MonitorIntegrationWavelengthMax", 0.0);
    TS_ASSERT_THROWS(alg->execute(), std::invalid_argument);
  }

  void test_execute_one_tranmission_1() {
    // One transmission run
    // Monitors : Yes
    // Monitor background : Yes
    // Monitor integration : Yes
    // Processing instructions : 0

    IAlgorithm_sptr alg =
        AlgorithmManager::Instance().create("CreateTransmissionWorkspace");

    alg->setChild(true);
    alg->initialize();

    alg->setProperty("FirstTransmissionRun", m_pointDetectorWS);
    alg->setProperty("WavelengthMin", 1.0);
    alg->setProperty("WavelengthMax", 15.0);
    alg->setProperty("I0MonitorIndex", 1);
    alg->setProperty("MonitorBackgroundWavelengthMin", 1.0);
    alg->setProperty("MonitorBackgroundWavelengthMax", 2.0);
    alg->setProperty("MonitorIntegrationWavelengthMin", 4.0);
    alg->setProperty("MonitorIntegrationWavelengthMax", 10.0);
    alg->setPropertyValue("ProcessingInstructions", "0");
    alg->setPropertyValue("OutputWorkspace", "demo_ws");
    alg->execute();
    MatrixWorkspace_sptr outWS = alg->getProperty("OutputWorkspace");

    std::cout << "outWS->getNumberHistograms() = " << outWS->getNumberHistograms() << "\n";
    std::cout << "outWS->blocksize() = "           << outWS->blocksize()           << "\n";
    std::cout << "outWS->readX(0)[0] = "           << outWS->readX(0)[0]           << "\n";
    std::cout << "outWS->readX(0)[7] = "           << outWS->readX(0)[7]           << "\n";
    std::cout << "outWS->readY(0)[0] = "           << outWS->readY(0)[0]           << "\n";
    std::cout << "outWS->readY(0)[7] = "           << outWS->readY(0)[7]           << "\n";

    /*
    TS_ASSERT_EQUALS(outWS->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outWS->blocksize(), 8);
    TS_ASSERT(outWS->readX(0)[0] >= 1.5);
    TS_ASSERT(outWS->readX(0)[7] <= 15.0);
    TS_ASSERT_DELTA(outWS->readY(0)[0], 9.4590, 0.0001);
    TS_ASSERT_DELTA(outWS->readY(0)[7], 9.4590, 0.0001);*/
    TS_ASSERT_EQUALS("Wavelength", outWS->getAxis(0)->unit()->unitID());
  }

  void test_execute_one_tranmission_2() {
    // One transmission run
    // Monitors : No
    // Monitor background : No
    // Monitor integration : No
    // Analysis mode : multi detector

  }

  void test_execute_one_transmission_3() {
    // One transmission run
    // Monitors : Yes
    // Monitor background : No
    // Monitor integration : No
    // Analysis mode : multi detector

  }

  void test_execute_one_transmission_4() {
    // One transmission run
    // Monitors : Yes
    // Monitor background : Yes
    // Monitor integration : No
    // Analysis mode : multi detector

  }

  void test_execute_one_transmission_5() {
    // One transmission run
    // Monitors : Yes
    // Monitor background : Yes
    // Monitor integration : Yes
    // Analysis mode : multi detector

  }

  void test_execute_two_transmission_runs() {

    // TODO: Use two transmission runs and test some values
  }

  // TODO: add more tests if needed

};

#endif /* ALGORITHMS_TEST_CREATETRANSMISSIONWORKSPACETEST_H_ */
