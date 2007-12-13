#ifndef LOADINSTRUMENTTEST_H_
#define LOADINSTRUMENTTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidDataHandling/LoadInstrument.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/Instrument.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidAPI/AnalysisDataService.h"
#include "MantidKernel/Exception.h"
#include "MantidAPI/FrameworkManager.h"
#include "MantidAPI/Workspace.h"
#include "MantidAPI/Algorithm.h"
#include "Component.h"
#include <vector>

using namespace Mantid::API;
using namespace Mantid::Kernel;
using namespace Mantid::DataHandling;
using Mantid::DataObjects::Workspace2D;

class LoadInstrumentTest : public CxxTest::TestSuite
{
public: 
  
  LoadInstrumentTest()
  {	
	//initialise framework manager to allow logging
	FrameworkManager manager;
	manager.initialize();
  }
  void testInit()
  {
    TS_ASSERT( !loader.isInitialized() );
    TS_ASSERT_THROWS_NOTHING(loader.initialize());    
    TS_ASSERT( loader.isInitialized() );
  }
  
  void testExec()
  {
    if ( !loader.isInitialized() ) loader.initialize();

    //create a workspace with some sample data
    wsName = "LoadInstrumentTest";
    WorkspaceFactory *factory = WorkspaceFactory::Instance();
    Workspace *ws = factory->create("Workspace2D");
    Workspace2D *ws2D = dynamic_cast<Workspace2D*>(ws);
    int histogramNumber = 2584;
    int timechannels = 100;
    ws2D->setHistogramNumber(histogramNumber);
    //loop to create data
    for (int i = 0; i < 2584; i++)
    {
      std::vector<double> timeChannelsVec(timechannels);
      std::vector<double> v(timechannels);
      // Create and fill another vector for the errors
      std::vector<double> e(timechannels);    
    //timechannels
      for (int j = 0; j < timechannels; j++)
      {
        timeChannelsVec[j] = j*100;
        v[j] = (i+j)%256;   
        e[j] = (i+j)%78;
      }
      // Populate the workspace.
      ws2D->setX(i, timeChannelsVec);
      ws2D->setData(i, v, e);
    }

    //put this workspace in the data service
    AnalysisDataService *data = AnalysisDataService::Instance();
    TS_ASSERT_THROWS_NOTHING(data->add(wsName, ws2D));    

    // Path to test input file assumes Test directory checked out from SVN
    inputFile = "../../../../Test/Instrument/HET_Definition.txt";
    loader.setProperty("Filename", inputFile);

    loader.setProperty("Workspace", wsName);
    
    std::string result;
    TS_ASSERT_THROWS_NOTHING( result = loader.getPropertyValue("Filename") )
    TS_ASSERT( ! result.compare(inputFile));

    TS_ASSERT_THROWS_NOTHING( result = loader.getPropertyValue("Workspace") )
    TS_ASSERT( ! result.compare(wsName));


    TS_ASSERT_THROWS_NOTHING(loader.execute());    

    TS_ASSERT( loader.isExecuted() );    
    
    // Get back the saved workspace
    Workspace *output;
    TS_ASSERT_THROWS_NOTHING(output = data->retrieve(wsName));
    
    Instrument& i = output->getInstrument();
    Mantid::Geometry::Component* source = i.getSource();
    TS_ASSERT_EQUALS( source->getName(), "Source");
    TS_ASSERT_EQUALS( source->getPos(), Mantid::Geometry::V3D(0,0,0));

    Mantid::Geometry::Component* samplepos = i.getSamplePos();
    TS_ASSERT_EQUALS( samplepos->getName(), "SamplePos");
    TS_ASSERT_EQUALS( samplepos->getPos(), Mantid::Geometry::V3D(0,10,0));

    TS_ASSERT_EQUALS(i.getDetectors()->nelements(),2184);

    Mantid::Geometry::Detector *ptrDet1000 = i.getDetector(1000);
    TS_ASSERT_EQUALS( ptrDet1000->getID(), 1000);
    TS_ASSERT_EQUALS( ptrDet1000->getName(), "PSD");
    TS_ASSERT_DELTA( ptrDet1000->getPos().X(), 3.86,0.01);
    TS_ASSERT_DELTA( ptrDet1000->getPos().Y(), 11.12,0.01);
    TS_ASSERT_DELTA( ptrDet1000->getPos().Z(), 0.43,0.01);
    TS_ASSERT_EQUALS( ptrDet1000->type(), "DetectorComponent");

    // Test input data is unchanged
    Workspace2D *output2DInst = dynamic_cast<Workspace2D*>(output);
    // Should be 2584 
    TS_ASSERT_EQUALS( output2DInst->getHistogramNumber(), histogramNumber);

  }
  
  void testFinal()
  {
    if ( !loader.isInitialized() ) loader.initialize();
    
    // The final() method doesn't do anything at the moment, but test anyway
    TS_ASSERT_THROWS_NOTHING(loader.finalize());    
    TS_ASSERT( loader.isFinalized() );
  }
  
private:
  LoadInstrument loader;
  std::string inputFile;
  std::string wsName;
  
};
  
#endif /*LOADINSTRUMENTTEST_H_*/
