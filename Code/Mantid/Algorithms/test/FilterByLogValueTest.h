/*
 * FilterByLogValueTest.h
 *
 *  Created on: Sep 15, 2010
 *      Author: janik
 */

#ifndef FILTERBYLOGVALUETEST_H_
#define FILTERBYLOGVALUETEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidAlgorithms/FilterByLogValue.h"
#include "MantidDataHandling/LoadEventPreNeXus.h"
#include "MantidKernel/DateAndTime.h"
#include "MantidKernel/TimeSeriesProperty.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidAPI/AlgorithmManager.h"

using namespace Mantid::Algorithms;
using namespace Mantid::DataHandling;
using namespace Mantid::DataObjects;
using namespace Mantid::Kernel;
using namespace Mantid::API;

class FilterByLogValueTest : public CxxTest::TestSuite
{
public:
  FilterByLogValueTest()
  {
    inputWS = "eventWS";
  }


  /** Setup for loading raw data */
  void setUp_Event()
  {
    IAlgorithm_sptr loader = AlgorithmManager::Instance().create("LoadSNSEventNexus");
    loader->initialize();
    loader->setPropertyValue("Filename", "../../../../Test/AutoTestData/CNCS_7850_event.nxs");
    loader->setPropertyValue("OutputWorkspace", inputWS);
    loader->setPropertyValue("filterByTof_Min", "45000");
    loader->setPropertyValue("filterByTof_Min", "50000");
    loader->execute();
    TS_ASSERT (loader->isExecuted() );

//    LoadEventPreNeXus * loader;
//    loader.initialize();
//    std::string eventfile( "../../../../Test/AutoTestData/CNCS_12772_neutron_event.dat" );
//    std::string pulsefile( "../../../../Test/AutoTestData/CNCS_12772_pulseid.dat" );
//    loader.setPropertyValue("EventFilename", eventfile);
//    loader.setProperty("PulseidFilename", pulsefile);
//    loader.setPropertyValue("MappingFilename", "../../../../Test/AutoTestData/CNCS_TS_2008_08_18.dat");
//    loader.setPropertyValue("PadEmptyPixels", "0");
    //    loader.setPropertyValue("OutputWorkspace", inputWS);
//    loader.execute();
//    TS_ASSERT (loader.isExecuted() );
  }



  void doTest(std::string outputWS)
  {
    //Retrieve Workspace
    this->setUp_Event();
    WS = boost::dynamic_pointer_cast<EventWorkspace>(AnalysisDataService::Instance().retrieve(inputWS));
    TS_ASSERT( WS ); //workspace is loaded

    size_t start_blocksize = WS->blocksize();
    size_t num_events = WS->getNumberEvents();
    double start_proton_charge = WS->run().getProtonCharge();
    size_t num_sample_logs = WS->run().getProperties().size();
    TS_ASSERT_EQUALS( num_events, 1208875 );

    //Do the filtering now.
    FilterByLogValue * alg = new FilterByLogValue();
    alg->initialize();
    alg->setPropertyValue("InputWorkspace", inputWS);
    alg->setPropertyValue("OutputWorkspace", outputWS);
    alg->setPropertyValue("LogName", "proton_charge");
    //We set the minimum high enough to cut out some real charge too, not just zeros.
    alg->setPropertyValue("MinimumValue", "1.85e7");
    alg->setPropertyValue("MaximumValue", "1e20");
    alg->setPropertyValue("TimeTolerance", "4e-2");

    alg->execute();
    TS_ASSERT( alg->isExecuted() );

    //Retrieve Workspace changed
    EventWorkspace_sptr outWS;
    outWS = boost::dynamic_pointer_cast<EventWorkspace>(AnalysisDataService::Instance().retrieve(outputWS));
    TS_ASSERT( outWS ); //workspace is loaded

    //Things that haven't changed
    TS_ASSERT_EQUALS( outWS->blocksize(), WS->blocksize());
    TS_ASSERT_EQUALS( outWS->getNumberHistograms(), WS->getNumberHistograms());

    //There should be some events
    TS_ASSERT_LESS_THAN( 0, outWS->getNumberEvents());

    TS_ASSERT_LESS_THAN( outWS->getNumberEvents(), num_events);
    TS_ASSERT_DELTA(  outWS->getNumberEvents() , 6536, 100);

    //Proton charge is lower
    TS_ASSERT_EQUALS( outWS->run().getProperties().size(), num_sample_logs);
    TS_ASSERT_LESS_THAN( outWS->run().getProtonCharge(), start_proton_charge );
    // But not 0
    TS_ASSERT_LESS_THAN( 0, outWS->run().getProtonCharge());

    //Still has a spectraDetectorMap;
    outWS->spectraMap();

  }

  void test_exec_renamed()
  {
    doTest(inputWS + "_filtered");
    AnalysisDataService::Instance().remove(inputWS);
    AnalysisDataService::Instance().remove(inputWS + "_filtered");
  }

  void test_exec_inplace()
  {
    doTest(inputWS);
    AnalysisDataService::Instance().remove(inputWS);
  }


  void setUp_Event2()
  {
    inputWS = "eventWS2";
    LoadEventPreNeXus loader;
    loader.initialize();
    std::string eventfile( "../../../../Test/AutoTestData/CNCS_7850_neutron_event.dat" );
    std::string pulsefile( "../../../../Test/AutoTestData/CNCS_7850_pulseid.dat" );
    loader.setPropertyValue("EventFilename", eventfile);
    loader.setProperty("PulseidFilename", pulsefile);
    loader.setPropertyValue("MappingFilename", "../../../../Test/AutoTestData/CNCS_TS_2008_08_18.dat");
    loader.setPropertyValue("OutputWorkspace", inputWS);
    loader.execute();
    TS_ASSERT (loader.isExecuted() );
  }

  void xtestExec2_slow()
  {
    std::string outputWS;
    this->setUp_Event2();

    //Retrieve Workspace
    WS = boost::dynamic_pointer_cast<EventWorkspace>(AnalysisDataService::Instance().retrieve(inputWS));
    TS_ASSERT( WS ); //workspace is loaded
    size_t start_blocksize = WS->blocksize();
    size_t num_events = WS->getNumberEvents();

    //Do the filtering now.
    FilterByLogValue * alg = new FilterByLogValue();
    alg->initialize();
    alg->setPropertyValue("InputWorkspace", inputWS);
    outputWS = "eventWS_relative";
    alg->setPropertyValue("OutputWorkspace", outputWS);
    alg->setPropertyValue("LogName", "proton_charge");
    //We set the minimum high enough to cut out some real charge too, not just zeros.
    alg->setPropertyValue("MinimumValue", "5e6");
    alg->setPropertyValue("MaximumValue", "1e20");
    alg->setPropertyValue("TimeTolerance", "3e-3");

    alg->execute();
    TS_ASSERT( alg->isExecuted() );

    //Retrieve Workspace changed
    EventWorkspace_sptr outWS;
    outWS = boost::dynamic_pointer_cast<EventWorkspace>(AnalysisDataService::Instance().retrieve(outputWS));
    TS_ASSERT( outWS ); //workspace is loaded

    //Things that haven't changed
    TS_ASSERT_EQUALS( outWS->blocksize(), WS->blocksize());
    TS_ASSERT_EQUALS( outWS->getNumberHistograms(), WS->getNumberHistograms());

    //There should be some events
    TS_ASSERT_LESS_THAN( 0, outWS->getNumberEvents());

    // Not many events left: 34612
    TS_ASSERT_LESS_THAN( outWS->getNumberEvents(), WS->getNumberEvents());
    TS_ASSERT_DELTA(  outWS->getNumberEvents() , 1093284, 100);

    //Proton charge is lower
    TS_ASSERT_LESS_THAN( outWS->run().getProtonCharge(), WS->run().getProtonCharge() );

    //Check the log entries
    TimeSeriesProperty<double> * log = dynamic_cast<TimeSeriesProperty<double> * >(outWS->run().getProperty("proton_charge"));
    TS_ASSERT(log);
    for (std::size_t i=0; i<log->realSize(); i++)
      TS_ASSERT_LESS_THAN( 0, log->nthValue(i));

  }



private:
  std::string inputWS;
  EventWorkspace_sptr WS;


};


#endif


