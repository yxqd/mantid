/*WIKI* 

Loads the mapping table between spectra and [[IDetector]] from a RAW file. It fills the [[SpectraToDetectorMap]] object contained in a [[workspace]]. This algorithm will fail if the [[workspace]] does not already point to a full [[instrument]] [[geometry]] (which usually means it must be run after [[LoadInstrument]]/[[LoadInstrumentFromRaw]]).

The association is one to many, i.e. a spectrum can have one or many detectors contributing to it. Alternatively the same spectrum can contribute to different spectra (for example in DAE2 (Data Aquisition Electronic) when a spectra containing electronically focussed data is created simultaneously with individual spectra).


*WIKI*/
#include "MantidDataHandling/LoadMappingTable.h"
#include "LoadRaw/isisraw2.h"
#include "MantidAPI/FileProperty.h"
#include "MantidAPI/SpectrumDetectorMapping.h"

namespace Mantid
{
namespace DataHandling
{

using namespace Kernel;
using namespace API;

DECLARE_ALGORITHM(LoadMappingTable)

/// Sets documentation strings for this algorithm
void LoadMappingTable::initDocs()
{
  this->setWikiSummary("Builds up the mapping between spectrum number and the detector objects in the [[instrument]] [[Geometry]].");
  this->setOptionalMessage("Builds up the mapping between spectrum number and the detector objects in the instrument Geometry.");
}


LoadMappingTable::LoadMappingTable() : Algorithm()
{
}

void LoadMappingTable::init()
{
  declareProperty(new FileProperty("Filename","", FileProperty::Load),
    "The name of the RAW file from which to obtain the mapping information, including its full or relative path." );
  declareProperty(new WorkspaceProperty<>("Workspace","Anonymous",Direction::InOut),
   "The name of the input and output workspace on which to perform the algorithm.");
}

void LoadMappingTable::exec()
{
  //Get the raw file name
  m_filename = getPropertyValue("Filename");
  // Get the input workspace
  const MatrixWorkspace_sptr localWorkspace = getProperty("Workspace");
        
  /// ISISRAW class instance which does raw file reading. Shared pointer to prevent memory leak when an exception is thrown.
  boost::scoped_ptr<ISISRAW2> iraw(new ISISRAW2);

  if (iraw->readFromFile(m_filename.c_str(),0) != 0) // ReadFrom File with no data
  {
    g_log.error("Unable to open file " + m_filename);
    throw Kernel::Exception::FileError("Unable to open File:" , m_filename);
  }
  progress(0.5);
  const int number_spectra=iraw->i_det; // Number of entries in the spectra/udet table
  if ( number_spectra == 0 )
  {
    g_log.warning("The spectra to detector mapping table is empty");
  }
  // Fill in the mapping in the workspace's ISpectrum objects
  localWorkspace->updateSpectraUsing(SpectrumDetectorMapping(iraw->spec,iraw->udet,number_spectra));
  progress(1);

  return;
}

} // Namespace DataHandling
} // Namespace Mantid
