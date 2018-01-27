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


def StartLiveData(*args, **kwargs):
    """
    StartLiveData dynamically adds the properties of the specific LiveListener
    that is used to itself, to allow usage such as the following:

        StartLiveData(Instrument='ISIS_Histogram', ...
                      PeriodList=[1,3], SpectraList=[2,4,6])

    Where PeriodList and SpectraList are properties of the ISISHistoDataListener
    rather than of StartLiveData. For StartLiveData to know those are valid
    properties, however, it first needs to know what the Instrument is.

    This is a similar situation as in the Load algorithm, where the Filename
    must be provided before other properties become available, and so it is
    solved here in the same way.
    """
    instrument, = _get_mandatory_args('StartLiveData', ["Instrument"], *args, **kwargs)

    # Create and execute
    (_startProgress, _endProgress, kwargs) = extract_progress_kwargs(kwargs)
    algm = _create_algorithm_object('StartLiveData',
                                    startProgress=_startProgress,
                                    endProgress=_endProgress)
    _set_logging_option(algm, kwargs)
    _set_store_ads(algm, kwargs)

    # Some properties have side effects and must be set separately
    def handleSpecialProperty(name, value=None):
        try:
            if value is None:
                value = kwargs.pop(name)
            algm.setProperty(name, value)

        except ValueError as ve:
            raise ValueError('Problem when setting %s. This is the detailed error '
                             'description: %s' % (name, str(ve)))
        except KeyError:
            pass  # ignore if kwargs[name] doesn't exist

    # Listener properties depend on these values, so they must be set first
    handleSpecialProperty('Instrument', instrument)
    handleSpecialProperty('Connection')
    handleSpecialProperty('Listener')

    # LHS Handling currently unsupported for StartLiveData
    lhs = _kernel.funcinspect.lhs_info()
    if lhs[0] > 0:  # Number of terms on the lhs
        raise RuntimeError("Assigning the output of StartLiveData is currently "
                           "unsupported due to limitations of the simpleapi. "
                           "Please call StartLiveData without assigning it to "
                           "to anything.")

    lhs_args = _get_args_from_lhs(lhs, algm)
    final_keywords = _merge_keywords_with_lhs(kwargs, lhs_args)

    # Check for any properties that aren't known and warn they will not be used
    for key in list(final_keywords.keys()):
        if key not in algm:
            logger.warning("You've passed a property (%s) to StartLiveData() "
                           "that doesn't apply to this Instrument." % key)
            del final_keywords[key]

    set_properties(algm, **final_keywords)
    algm.execute()

    return _gather_returns("StartLiveData", lhs, algm)
