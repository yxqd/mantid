#  This file is part of the mantid workbench.
#
#  Copyright (C) 2018 mantidproject
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import (absolute_import, division, print_function, unicode_literals)




def Load(*args, **kwargs):
    """
    Load is a more flexible algorithm than other Mantid algorithms.
    It's aim is to discover the correct loading algorithm for a
    given file. This flexibility comes at the expense of knowing the
    properties out right before the file is specified.

    The argument list for the Load function has to be more flexible to
    allow this searching to occur. Two arguments must be specified:

      - Filename :: The name of the file,
      - OutputWorkspace :: The name of the workspace,

    either as the first two arguments in the list or as keywords. Any other
    properties that the Load algorithm has can be specified by keyword only.

    Some common keywords are:
     - SpectrumMin,
     - SpectrumMax,
     - SpectrumList,
     - EntryNumber

    Example:
      # Simple usage, ISIS NeXus file
      run_ws = Load('INSTR00001000.nxs')

      # Histogram NeXus with SpectrumMin and SpectrumMax = 1
      run_ws = Load('INSTR00001000.nxs', SpectrumMin=1,SpectrumMax=1)

      # Event NeXus with precount on
      event_ws = Load('INSTR_1000_event.nxs', Precount=True)

      # The output workspace name is picked up from the LHS unless overridden
      Load('INSTR00001000.nxs',OutputWorkspace='run_ws')
    """
    filename, = _get_mandatory_args('Load', ["Filename"], *args, **kwargs)
    if not filename:
        # If we try to set property with a None type we get a unhelpful error about allocators
        # so check up front here
        raise ValueError("Problem with supplied Filename. The value given was a 'None' "
                         "type and cannot be used. Please ensure the Filename is set"
                         " to the path of the file.")

    # Create and execute
    (_startProgress, _endProgress, kwargs) = extract_progress_kwargs(kwargs)
    algm = _create_algorithm_object('Load', startProgress=_startProgress,
                                    endProgress=_endProgress)
    _set_logging_option(algm, kwargs)
    _set_store_ads(algm, kwargs)
    try:
        algm.setProperty('Filename', filename)  # Must be set first
    except ValueError as ve:
        raise ValueError('Problem when setting Filename. This is the detailed error '
                         'description: ' + str(ve) + '\nIf the file has been found '
                         'but you got this error, you might not have read permissions '
                         'or the file might be corrupted.\nIf the file has not been found, '
                         'you might have forgotten to add its location in the data search '
                         'directories.')
    # Remove from keywords so it is not set twice
    if 'Filename' in kwargs:
        del kwargs['Filename']
    lhs = _kernel.funcinspect.lhs_info()
    # If the output has not been assigned to anything, i.e. lhs[0] = 0 and kwargs does not have OutputWorkspace
    # then raise a more helpful error than what we would get from an algorithm
    if lhs[0] == 0 and 'OutputWorkspace' not in kwargs:
        raise RuntimeError("Unable to set output workspace name. Please either assign the output of "
                           "Load to a variable or use the OutputWorkspace keyword.")

    lhs_args = _get_args_from_lhs(lhs, algm)
    final_keywords = _merge_keywords_with_lhs(kwargs, lhs_args)
    # Check for any properties that aren't known and warn they will not be used
    for key in list(final_keywords.keys()):
        if key not in algm:
            logger.warning("You've passed a property (%s) to Load() that doesn't apply to this file type." % key)
            del final_keywords[key]
    set_properties(algm, **final_keywords)
    algm.execute()

    # If a WorkspaceGroup was loaded then there will be a set of properties that have an underscore in the name
    # and users will simply expect the groups to be returned NOT the groups + workspaces.
    return _gather_returns('Load', lhs, algm, ignore_regex=['LoaderName', 'LoaderVersion', '.*_.*'])
