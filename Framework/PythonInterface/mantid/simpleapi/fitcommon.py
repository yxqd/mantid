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


def fitting_algorithm(inout=False):
    """
    Decorator generating code for fitting algorithms (Fit, CalculateChiSquared,
    EvaluateFunction).
    When applied to a function definition this decorator replaces its code
    with code of function 'wrapper' defined below.
    :param inout: if True, return also the InOut properties of algorithm f
    """
    def inner_fitting_algorithm(f):
        """
        :param f: algorithm calling Fit
        """
        def wrapper(*args, **kwargs):
            function, input_workspace = _get_mandatory_args(function_name,
                                                            ["Function", "InputWorkspace"],
                                                            *args, **kwargs)
            # Remove from keywords so it is not set twice
            if "Function" in kwargs:
                del kwargs['Function']
            if "InputWorkspace" in kwargs:
                del kwargs['InputWorkspace']

            # Check for behaviour consistent with old API
            if type(function) == str and function in _api.AnalysisDataService:
                msg = "Fit API has changed. The function must now come " + \
                      "first in the argument list and the workspace second."
                raise ValueError(msg)
            # Deal with case where function is a FunctionWrapper.
            if isinstance(function, FunctionWrapper):
                function = function.__str__()

            # Create and execute
            algm = _create_algorithm_object(function_name)
            _set_logging_option(algm, kwargs)
            _set_store_ads(algm, kwargs)
            if 'EvaluationType' in kwargs:
                algm.setProperty('EvaluationType', kwargs['EvaluationType'])
                del kwargs['EvaluationType']
            algm.setProperty('Function', function)  # Must be set first
            if input_workspace is not None:
                algm.setProperty('InputWorkspace', input_workspace)
            else:
                del algm['InputWorkspace']

            # Set all workspace properties before others
            for key in list(kwargs.keys()):
                if key.startswith('InputWorkspace_'):
                    algm.setProperty(key, kwargs[key])
                    del kwargs[key]

            lhs = _lhs_info()
            # Check for unknown properties and warn they will not be used
            for key in list(kwargs.keys()):
                if key not in algm:
                    msg = 'Property {} to {} does not apply to any of the ' +\
                          ' input workspaces'.format(key, function_name)
                    logger.warning(msg)
                    del kwargs[key]
            set_properties(algm, **kwargs)
            algm.execute()
            return _gather_returns(function_name, lhs, algm, inout=inout)
        # end
        function_name = f.__name__
        signature = ("\bFunction, InputWorkspace", "**kwargs")
        fwrapper = _customise_func(wrapper, function_name, signature, f.__doc__)
        if function_name not in __SPECIALIZED_FUNCTIONS__:
            __SPECIALIZED_FUNCTIONS__.append(function_name)
        return fwrapper
    return inner_fitting_algorithm


