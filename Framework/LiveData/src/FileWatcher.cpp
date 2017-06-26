#include "MantidLiveData/FileWatcher.h"
#include <iostream>

namespace Mantid {
namespace LiveData {
	/**
	* The constructor
	*/
	FileWatcher::FileWatcher(std::string path)
		: m_path(path), m_changed(false), dw(m_path, Poco::DirectoryWatcher::DW_FILTER_ENABLE_ALL, 2) {}

	std::string FileWatcher::get_path() {
		return m_path;
	}
	
	/**
	* Returns whether any files in the watched directory have changed.
	*/
	bool FileWatcher::hasChanged() {
		return m_changed;
	}

	/**
	* Read and empty the list of changed files .
	*/
	std::vector<Poco::File> FileWatcher::readChanged() {
		std::vector<Poco::File> changedFiles = m_changedFiles;
		m_changedFiles.clear();
		m_changed = false;
		return changedFiles;
	}

	/**
	* Start watching the directory.
	*/
	void FileWatcher::start() {
		dw.itemAdded += Poco::delegate(this, &FileWatcher::onItemAdded);
		dw.itemRemoved += Poco::delegate(this, &FileWatcher::onItemRemoved);
		dw.itemModified += Poco::delegate(this, &FileWatcher::onItemModified);
		dw.itemMovedFrom += Poco::delegate(this, &FileWatcher::onItemMovedFrom);
		dw.itemMovedTo += Poco::delegate(this, &FileWatcher::onItemMovedTo);
	}

	/**
	* The actions to perform when a file has been added.
	*/
	void FileWatcher::onItemAdded(const Poco::DirectoryWatcher::DirectoryEvent& ev)
	{
		m_changedFiles.push_back(ev.item);
		std::cout << "ADDED: item at " << ev.item.path() << std::endl;
		m_changed = true;
	}

	/**
	* The actions to perform when a file has been removed.
	*/
	void FileWatcher::onItemRemoved(const Poco::DirectoryWatcher::DirectoryEvent& ev)
	{
		std::cout << "REMOVED: item at " << ev.item.path() << std::endl;
		m_changed = true;
	}

	/**
	* The actions to perform when a file has been modified.
	*/
	void FileWatcher::onItemModified(const Poco::DirectoryWatcher::DirectoryEvent& ev)
	{
		std::cout << "MODIFIED: item at " << ev.item.path() << std::endl;
		m_changed = true;
	}

	// Never happens (on win)?
	void FileWatcher::onItemMovedFrom(const Poco::DirectoryWatcher::DirectoryEvent& ev)
	{
		std::cout << "MOVED FROM: item at " << ev.item.path() << std::endl;
		m_changed = true;
	}

	// Never happens (on win)?
	void FileWatcher::onItemMovedTo(const Poco::DirectoryWatcher::DirectoryEvent& ev)
	{
		std::cout << "MOVED TO: item at " << ev.item.path() << std::endl;
		m_changed = true;
	}

	// TODO set changed method
	
	// TODO read 

	// TODO white/blacklist files? (regex?)

} // namespace LiveData
} // namespace Mantid
