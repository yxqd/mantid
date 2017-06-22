#ifndef MANTID_LIVEDATA_FILEWATCHERTEST_H_
#define MANTID_LIVEDATA_FILEWATCHERTEST_H_

#include <cxxtest/TestSuite.h>
#include <fstream>
#include <stdio.h>

#include "MantidLiveData/FileWatcher.h"

static const std::string folderPath = "C:\\Instrument\\FileWatcherTestBed\\";
static std::vector<std::string> testFiles;
using Mantid::LiveData::FileWatcher;

class FileWatcherTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static FileWatcherTest *createSuite() { return new FileWatcherTest(); }
  static void destroySuite( FileWatcherTest *suite ) { delete suite; }

  void test_pathIsSet()
  {
	FileWatcher fw(folderPath);
    TS_ASSERT_EQUALS(fw.get_path(), folderPath);
  }
  
  void test_flagIsFalse(){
	  FileWatcher fw(folderPath);
	  TS_ASSERT_EQUALS(fw.hasChanged(), false);
  }

  void test_flagSetToTrueIfFileAdded() {
	  FileWatcher fw(folderPath);
	  fw.start();

	  std::ofstream testFile = openFile("example.txt");
	  testFile.close();

	  Sleep(1000);
	  TS_ASSERT_EQUALS(fw.hasChanged(), true);
  }
  std::ofstream openFile(std::string fileName) {
	  std::ofstream file;

	  // check doesnt exist
	  std::string filePath = folderPath + fileName;
	  testFiles.push_back(filePath);

	  file.open(filePath);
	  return file;
  }

  void tearDown() override {
	  // Remove all files created
	  for (auto const& filePath : testFiles) {
		  remove(filePath.c_str());
	  }
  }

};

#endif /* MANTID_LIVEDATA_FILEWATCHERTEST_H_ */