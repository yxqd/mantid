#ifndef MANTID_ALGORITHMS_GENERATEEVENTSFILTER_H_
#define MANTID_ALGORITHMS_GENERATEEVENTSFILTER_H_

#include "MantidKernel/System.h"
#include "MantidAPI/Algorithm.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidKernel/TimeSeriesProperty.h"
#include "MantidDataObjects/SplittersWorkspace.h"
#include "MantidAPI/ISplittersWorkspace.h"
#include "MantidAPI/ITableWorkspace.h"

namespace Mantid
{
namespace Algorithms
{

  /** GenerateEventsFilter : Generate an events-filter, i.e., a SplittersWorkspace according
      to user's request.

      Request can be a combination of log value and time, i.e.,

      (1) T_min
      (2) T_max
      (3) delta Time
      (4) Min log value
      (5) Max log value
      (6) delta log value
      (7) identify log value increment (bool)
      (8) number of sections per interval (applied to log value only!)

      This algorithm can generate filters including
      (1) deltaT per interval from T_min to T_max with : Log value is not given
      (2) delta log value per interval from T_min to T_max and from min log value to max log value

      Note:
      (1) Time can be (a) relative time in ns  (b) relative time in second (float) (c) percentage time
      (2) if option "identify log value increment"
    
    @date 2012-04-09

    Copyright &copy; 2012 ISIS Rutherford Appleton Laboratory & NScD Oak Ridge National Laboratory

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
  class DLLExport GenerateEventsFilter : public API::Algorithm
  {
  public:
    explicit GenerateEventsFilter();
    virtual ~GenerateEventsFilter();
    
    /// Algorithm's name for identification overriding a virtual method
    virtual const std::string name() const { return "GenerateEventsFilter";}
    /// Algorithm's version for identification overriding a virtual method
    virtual int version() const { return 1;}
    /// Algorithm's category for identification overriding a virtual method
    virtual const std::string category() const { return "Events\\EventFiltering";}

  private:
    /// Sets documentation strings for this algorithm
    virtual void initDocs();
    // Implement abstract Algorithm methods
    void init();
    // Implement abstract Algorithm methods
    void exec();

    /// Process properties
    void processInOutWorkspaces();

    void processInputTime();
    void setFilterByTimeOnly();
    void setFilterByLogValue(std::string logname);

    void processSingleValueFilter(double minvalue, double maxvalue,
        bool filterincrease, bool filterdecrease);

    void processMultipleValueFilters(double minvalue, double maxvalue,
        bool filterincrease, bool filterdecrease);

    void makeFilterByValue(Kernel::TimeSplitterType& split, double min, double max, double TimeTolerance, bool centre,
        bool filterIncrease, bool filterDecrease, Kernel::DateAndTime startTime, Kernel::DateAndTime stopTime,
        int wsindex);

    /// Make multiple-log-value filters in serial
    void makeMultipleFiltersByValues(std::map<size_t, int> indexwsindexmap, std::vector<double> logvalueranges, bool centre,
                                     bool filterIncrease, bool filterDecrease, Kernel::DateAndTime startTime,
                                     Kernel::DateAndTime stopTime);

    /// Make multiple-log-value filters in serial in parallel
    void makeMultipleFiltersByValuesParallel(std::map<size_t, int> indexwsindexmap, std::vector<double> logvalueranges, bool centre,
                                             bool filterIncrease, bool filterDecrease, Kernel::DateAndTime startTime,
                                             Kernel::DateAndTime stopTime);

    void makeMultipleFiltersByValuesPartialLog(int istart, int iend,
                                               std::vector<Kernel::DateAndTime>& vecSplitTime,
                                               std::vector<int>& vecSplitGroup,
                                               std::map<size_t, int> indexwsindexmap,
                                               const std::vector<double>& logvalueranges, Kernel::time_duration tol,
                                               bool filterIncrease, bool filterDecrease,
                                               Kernel::DateAndTime startTime, Kernel::DateAndTime stopTime);


    void processIntegerValueFilter(int minvalue, int maxvalue,
                                   bool filterIncrease, bool filterDecrease, Kernel::DateAndTime runend);

    size_t searchValue(std::vector<double> sorteddata, double value);

    /// Add a splitter
    void addSplitter(Kernel::DateAndTime starttime, Kernel::DateAndTime stoptime, int wsindex, std::string info);

    /// Add a splitter
    void make_splitter(Kernel::DateAndTime start, Kernel::DateAndTime stop, int group, Kernel::time_duration tolerance);

    /// Create a splitter and add to the vector of time splitters
    Kernel::DateAndTime makeSplitterInVector(std::vector<Kernel::DateAndTime>& vecSplitTime, std::vector<int>& vecGroupIndex,
                                             Kernel::DateAndTime start, Kernel::DateAndTime stop, int group,
                                             Kernel::time_duration tolerance, int64_t tol_ns, Kernel::DateAndTime lasttime);


    /// Generate a matrix workspace containing splitters
    void generateSplittersInMatrixWorkspace();

    void generateSplittersInMatrixWorkspaceParallel();

    DataObjects::EventWorkspace_const_sptr m_dataWS;

    /// SplitterWorkspace
    API::ISplittersWorkspace_sptr m_splitWS;
    /// Matrix workspace containing splitters
    API::MatrixWorkspace_sptr m_filterWS;

    API::ITableWorkspace_sptr m_filterInfoWS;

    Kernel::DateAndTime m_startTime;
    Kernel::DateAndTime m_stopTime;

    double m_timeUnitConvertFactorToNS;

    Kernel::TimeSeriesProperty<double>* m_dblLog;
    Kernel::TimeSeriesProperty<int>* m_intLog;

    bool m_logAtCentre;
    double m_logTimeTolerance;

    /// Flag to output matrix workspace for fast log
    bool m_forFastLog;

    /// SplitterType
    Kernel::TimeSplitterType m_splitters;
    /// Vector as date and time
    std::vector<Kernel::DateAndTime> m_vecSplitterTime;
    std::vector<int> m_vecSplitterGroup;

    /// Processing algorithm type
    bool m_useParallel;

    std::vector<std::vector<Kernel::DateAndTime> > vecSplitterTimeSet;
    std::vector<std::vector<int> > vecGroupIndexSet;

  };

  /** Generate a new time splitter and add to a list of splitters

  void make_splitter(Kernel::DateAndTime start, Kernel::DateAndTime stop, int group, Kernel::time_duration tolerance,
                     Kernel::TimeSplitterType& splitters)
  {
    Kernel::SplittingInterval newsplit(start - tolerance, stop - tolerance, group);
    splitters.push_back(newsplit);

    return;
  }
  */

} // namespace Algorithms
} // namespace Mantid

#endif  /* MANTID_ALGORITHMS_GENERATEEVENTSFILTER_H_ */
