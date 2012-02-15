#ifndef MANTID_VATES_MD_LOADING_PRESENTER
#define MANTID_VATES_MD_LOADING_PRESENTER

#include "MantidKernel/System.h"
#include <vtkDataSet.h>
#include <string>
#include <vector>

class vtkUnstructuredGrid;
namespace Mantid
{
  namespace VATES
  {
    class ProgressAction;
    class vtkDataSetFactory;
     /** 
    @class MDLoadingPresenter
    Abstract presenters for loading conversion of MDEW workspaces into render-able vtk objects.
    @author Owen Arnold, Tessella plc
    @date 05/08/2011

    Copyright &copy; 2011 ISIS Rutherford Appleton Laboratory & NScD Oak Ridge National Laboratory

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

    File change history is stored at: <https://svn.mantidproject.org/mantid/trunk/Code/Mantid>
    Code Documentation is available at: <http://doxygen.mantidproject.org>
    */
    class DLLExport MDLoadingPresenter
    {
    public:
      virtual vtkDataSet* execute(vtkDataSetFactory* factory, ProgressAction& rebinningProgressUpdate, ProgressAction& drawingProgressUpdate) = 0;
      virtual void executeLoadMetadata() = 0;
      virtual bool hasTDimensionAvailable() const = 0;
      virtual std::vector<double> getTimeStepValues() const = 0;
      virtual bool canReadFile() const = 0;
      virtual const std::string& getGeometryXML() const = 0;
      virtual ~MDLoadingPresenter(){}
      virtual std::string getWorkspaceTypeName()
      {
        return "NotSet";
      }
    };

  }
}

#endif
