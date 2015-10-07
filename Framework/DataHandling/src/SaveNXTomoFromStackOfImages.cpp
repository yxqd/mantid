#include "MantidDataHandling/SaveNXTomoFromStackOfImages.h"

namespace Mantid {
namespace DataHandling {

using Mantid::Kernel::Direction;
using Mantid::API::WorkspaceProperty;

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(SaveNXTomoFromStackOfImages)

/// Default constructor
SaveNXTomoFromStackOfImages::SaveNXTomoFromStackOfImages() {}

/// Destructor
SaveNXTomoFromStackOfImages::~SaveNXTomoFromStackOfImages() {}

/// Algorithms name for identification. @see Algorithm::name
const std::string SaveNXTomoFromStackOfImages::name() const {
  return "SaveNXTomoFromStackOfImages";
}

/// Algorithm's version for identification. @see Algorithm::version
int SaveNXTomoFromStackOfImages::version() const { return 1; }

/// Algorithm's category for identification. @see Algorithm::category
const std::string SaveNXTomoFromStackOfImages::category() const {
  return "DataHandling\\Nexus;DataHandling\\Tomography;Diffraction";
}

/// Algorithm's summary for use in the GUI and help. @see Algorithm::summary
const std::string SaveNXTomoFromStackOfImages::summary() const {
  return "Produces an NXTomo file from a set (stack) of individual image "
         "files.";
}

/**
 * Initialize the algorithm's properties.
 */
void SaveNXTomoFromStackOfImages::init() {

  declareProperty(new API::FileProperty("SampleFilename", "",
                                        API::FileProperty::Load, exts),
                  "First sample file in the stack of images/files. There must "
                  "be at least one.");

  auto mustBePositive = boost::make_shared<BoundedValidator<int64_t>>();
  declareProperty("SampleMaxIndex", EMPTY_INT(), mustBePositive,
                  "Maximum sample image index to load.");

  declareProperty(
      new API::FileProperty("FlatFieldFilename", "",
                            API::FileProperty::OptionalLoad, exts),
      "First flat or open beam image file. There can be none, one or many.");

  declareProperty(
      "FlatFieldMaxIndex", EMPTY_INT(), mustBePositive,
      "Maximum index to load when loading multiple flat field images.");

  declareProperty(
      new API::FileProperty("DarkFieldFilename", "",
                            API::FileProperty::OptionalLoad, exts),
      "First dark or open beam image file. There can be none, one or many.");

  declareProperty(
      "DarkFieldMaxIndex", EMPTY_INT(), mustBePositive,
      "Maximum index to load when loading multiple dark field images.");

  declareProperty(new PropertyWithValue<bool>("OverwriteFile", false,
                                              Kernel::Direction::Input),
                  "Replace any existing file of the same name instead of "
                  "appending data. If the file exists and this option is left "
                  "or given as false the algorithm will not do anything");

  declareProperty(new API::FileProperty("OutputFilename", "",
                                        API::FileProperty::Load, exts),
                  "Name of the output NXtomo. The extension .nxs will be "
                  "appended if missing.");
}

/**
 * Execute the algorithm.
 */
void SaveNXTomoFromStackOfImages::exec() {
  // TODO
}

} // namespace DataHandling
} // namespace Mantid
