#ifndef MANTID_LIVEDATA_FILEWATCHERTEST_H_
#define MANTID_LIVEDATA_FILEWATCHERTEST_H_

#include <cxxtest/TestSuite.h>
#include <fstream>
#include <stdio.h>

#include "MantidLiveData/FileWatcher.h"

static const std::string folderPath = "C:\\Instrument\\FileWatcherTestBed\\";
static std::vector<std::string> testFiles;
using Mantid::LiveData::FileWatcher;
using Poco::DirectoryWatcher;

static bool flag;

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

  void test_added()
  {
	  // Create filewatcher
	  std::string filePath = folderPath + "test_added.txt";
	  FileWatcher fw(folderPath);
	  fw.start();

	  Poco::Thread::sleep(1000);

	  // Add file in watched directory
	  Poco::FileOutputStream fos(filePath);
	  testFiles.push_back(filePath);
	  fos.close();

	  Poco::Thread::sleep(1000);

	  // Check updated flag has changed
	  TS_ASSERT_EQUALS(true, fw.hasChanged());

	  // Read changed files
	  std::set<Poco::File> addedFiles = fw.readChanges().added;
	  Poco::File added = *addedFiles.begin();
	  TS_ASSERT_EQUALS(filePath, added.path());
	  TS_ASSERT_EQUALS(false, fw.hasChanged());
  }

  void tearDown() override {
	  // Remove all files created
	  for (auto const& filePath : testFiles) {
		  remove(filePath.c_str());
	  }
	  testFiles.clear();
  }

};

#endif /* MANTID_LIVEDATA_FILEWATCHERTEST_H_ */