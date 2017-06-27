#ifndef MANTID_LIVEDATA_FILEWATCHER_H_
#define MANTID_LIVEDATA_FILEWATCHER_H_

#include "Poco/DirectoryWatcher.h"
#include "Poco/Delegate.h"
#include "Poco/FileStream.h"
#include "MantidKernel/System.h"

namespace Mantid {
namespace LiveData {

/** FileWatcher : TODO: DESCRIPTION

  Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
  National Laboratory & European Spallation Source

  This file is part of Mantid.

  Mantid is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  Mantid is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  File change history is stored at: <https://github.com/mantidproject/mantid>
  Code Documentation is available at: <http://doxygen.mantidproject.org>
*/
class DLLExport FileWatcher {

public:
	struct ChangedFiles {
		std::set<Poco::File> added;
		std::set<Poco::File> modified;
		std::set<Poco::File> deleted;

		void clear() {
			added.clear();
			modified.clear();
			deleted.clear();
		}
	};

	FileWatcher(std::string path);
	//~FileWatcher();
	std::string get_path();
	bool hasChanged();
	FileWatcher::ChangedFiles readChanges();

	void start();

private:
	bool m_changed;
	std::string m_path;
	Poco::DirectoryWatcher dw;
	struct ChangedFiles m_changedFiles;

	void onItemAdded(const Poco::DirectoryWatcher::DirectoryEvent& ev);
	void onItemRemoved(const Poco::DirectoryWatcher::DirectoryEvent & ev);
	void onItemModified(const Poco::DirectoryWatcher::DirectoryEvent & ev);
	void onItemMovedFrom(const Poco::DirectoryWatcher::DirectoryEvent & ev);
	void onItemMovedTo(const Poco::DirectoryWatcher::DirectoryEvent & ev);
};

} // namespace LiveData
} // namespace Mantid

#endif /* MANTID_LIVEDATA_FILEWATCHER_H_ */