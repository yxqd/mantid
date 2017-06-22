#include "MantidLiveData/FileWatcher.h"
#include <iostream>

namespace Mantid {
namespace LiveData {
	/**
	* The constructor
	*/
	FileWatcher::FileWatcher(std::string path)
		: m_path(path), m_changed(false) {}

	std::string FileWatcher::get_path() {
		return m_path;
	}

	bool FileWatcher::hasChanged() {
		return m_changed;
	}

	void FileWatcher::start() {
		Poco::DirectoryWatcher dw(m_path, Poco::DirectoryWatcher::DW_FILTER_ENABLE_ALL, 2);
		dw.itemAdded += Poco::delegate(this, &FileWatcher::onItemAdded);
		dw.itemRemoved += Poco::delegate(this, &FileWatcher::onItemRemoved);
		dw.itemModified += Poco::delegate(this, &FileWatcher::onItemModified);
		dw.itemMovedFrom += Poco::delegate(this, &FileWatcher::onItemMovedFrom);
		dw.itemMovedTo += Poco::delegate(this, &FileWatcher::onItemMovedTo);
		while (true) {}
	}

	void FileWatcher::onItemAdded(const Poco::DirectoryWatcher::DirectoryEvent& ev)
	{
		std::cout << "ADDED: item at " << ev.item.path() << std::endl;
		m_changed = true;
	}

	void FileWatcher::onItemRemoved(const Poco::DirectoryWatcher::DirectoryEvent& ev)
	{
		std::cout << "REMOVED: item at " << ev.item.path() << std::endl;
		m_changed = true;
	}

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
