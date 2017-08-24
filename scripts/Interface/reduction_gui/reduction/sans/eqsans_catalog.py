#pylint: disable=invalid-name
"""
    Data catalog for EQSANS
"""
from __future__ import (absolute_import, division, print_function)
from reduction_gui.reduction.sans.data_cat import DataCatalog as BaseCatalog
from reduction_gui.reduction.sans.data_cat import DataSet,  DataType
import re
import datetime
import traceback
from twisted.python.log import logerr

# Check whether Mantid is available
try:
    from mantid.api import AnalysisDataService
    import mantid.simpleapi as api
    HAS_MANTID = True
except:
    HAS_MANTID = False

try:
    import mantidplot
    from mantid.kernel import logger
    IN_MANTIDPLOT = True
except:
    IN_MANTIDPLOT = False
    import logging
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger("data_cat")


class EQSANSDataType(DataType):
    TABLE_NAME="eqsans_datatype"


class EQSANSDataSet(DataSet):
    TABLE_NAME="eqsans_dataset"
    data_type_cls = EQSANSDataType

    def __init__(self, run_number, title, run_start, duration, sdd):
        super(EQSANSDataSet, self).__init__(run_number, title, run_start, duration, sdd)

    
    @classmethod
    def handle(cls, file_path):
        """
            Return a DB handle for the given file, such as a run number
            This will handle file formats in two formats:
            EQSANS_([0-9]+)_event
            EQSANS_([0-9]+).nxs
        """
        file_path = file_path.strip()
        r_re = re.search("EQSANS_([0-9]+)(_event|\.nxs)", file_path)
        if r_re is not None:
            return r_re.group(1)
        else:
            # Check whether we simply have a run number
            try:
                int(file_path)
                return file_path
            except:
                return None
        return None
    
    
    @classmethod
    def load_meta_data(cls, file_path, run, out_ws_name):
        
        def read_prop(ws_object, prop):
            try:
                return str(ws_object.getRun().getProperty(prop).value)
            except:
                return ""

        def read_series(ws_object, prop):
            try:
                return float(ws_object.getRun().getProperty(prop).getStatistics().mean)
            except:
                return -1
        
        file_path = str(file_path)
        out_ws_name = str(out_ws_name)
        
        api.LoadEventNexus(
            Filename=file_path,
            OutputWorkspace=out_ws_name,
            MetaDataOnly=True
        )
        out_ws = api.mtd[out_ws_name]
        
        
        runno = read_prop(out_ws, "run_number")
        if runno=="":
            runno = run

        title = read_prop(out_ws, "run_title")
        t_str = read_prop(out_ws, "start_time")
        # Get rid of the training microseconds
        toks = t_str.split('.')
        if len(toks)>=2:
            t_str=toks[0]
        t = datetime.datetime.strptime(t_str, '%Y-%m-%dT%H:%M:%S')
        # TZ offset
        offset = datetime.datetime.now()-datetime.datetime.utcnow()
        t = t+offset
        run_start = t.strftime('%y-%m-%d %H:%M')

        duration = read_prop(out_ws, "duration")
        try:
            duration = float(duration)
        except:
            duration = 0

        sdd = read_series(out_ws, "detectorZ")

        api.DeleteWorkspace(out_ws_name)
        
        d = EQSANSDataSet(runno, title, run_start, duration, sdd)
        return d


class DataCatalog(BaseCatalog):
    extension = ["nxs", "nxs.h5"]
    data_set_cls = EQSANSDataSet

    def __init__(self, replace_db=False):
        super(DataCatalog, self).__init__(replace_db)
