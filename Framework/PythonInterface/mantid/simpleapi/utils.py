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
"""General utilities for creating and executing algorithms
"""
from __future__ import (absolute_import, division, print_function, unicode_literals)

# std imports
from collections import namedtuple, OrderedDict
import inspect

# 3rd party imports
from six import iteritems

# local imports
from mantid.api import (
    AlgorithmManagerImpl as AlgorithmManager,
    AnalysisDataServiceImpl as AnalysisDataService,
    FunctionProperty, IWorkspaceProperty
)
from mantid.kernel import DataItem, Direction, Property
from .functionwrapper import FunctionWrapper


PARENT_ALGORITHM_FN = "PyExec"


def create_algorithm(name, version, progress_start=0.0, progress_end=1.0,
                     enable_logging=True, store_in_ads=True):
    """
    Create and initialize the named algorithm of the given version. This
    method checks whether the function call has come from within a PyExec
    call. If that is the case then an unmanaged child algorithm is created.

    :param name: Name of the algorithm
    :param version: Version of algorithm required
    :param progress_start: An optional position to start progress reporting
    :param progress_end: An optional position to end progress reporting
    :param enable_logging: If true then enable logging for the algorithm
    :param store_in_ads: If true store the resultant workspaces in the ADS
    :return: A new instance of the named algorithm
    """
    parent = find_parent_algorithm(inspect.currentframe())
    if parent is not None:
        kwargs = dict(version=version)
        kwargs['startProgress'] = float(progress_start)
        kwargs['endProgress'] = float(progress_end)
        alg = parent.createChildAlgorithm(name, **kwargs)
        alg.setLogging(parent.isLogging())  # default is to log if parent is logging
    else:
        # managed algorithm so that progress reporting
        # can be more easily wired up automatically
        alg = AlgorithmManager.Instance().create(name, version)

    # set attributes
    alg.setRethrows(True)
    alg.setLogging(enable_logging)
    alg.setAlwaysStoreInADS(store_in_ads)

    return alg


def find_parent_algorithm(frame):
    """
    Look for a PyExec method in the call stack and return
    the self object that the method is attached to

    :param frame The starting frame for the stack walk
    :returns The self object that is running the PyExec method
             or None if one was not found
    """
    # Return the 'self' object of a given frame
    def get_self(frame_arg):
        return frame_arg.f_locals['self']

    # Look recursively for the method in the stack
    if frame.f_code.co_name == PARENT_ALGORITHM_FN:
        return get_self(frame)
    while True:
        if frame.f_back:
            if frame.f_back.f_code.co_name == PARENT_ALGORITHM_FN:
                return get_self(frame.f_back)
            frame = frame.f_back
        else:
            break
    if frame.f_code.co_name == PARENT_ALGORITHM_FN:
        return get_self(frame)
    else:
        return None


def merge_lhs_args(algm, lhs, kwargs):
    """
    Merge identifiers provided on the LHS with other keyword arguments
    to provide a final keyword dict for execution

    :param algm: An algorithm reference
    :param lhs: A list of variables given on the lhs
    :param kwargs: An existing set of keyword arguments provided in a function call
    :return: An updated keyword dict
    """
    lhs_len = len(lhs)
    if lhs_len == 0:
        return kwargs

    output_props = [algm.getProperty(p) for p in algm.outputProperties()]
    nprops = len(output_props)
    lhs_kwargs = {}
    lhs_idx = 0
    for p in output_props:
        if is_workspace_property(p):
            if 0 < lhs_len < nprops:
                lhs_kwargs[p.name] = ret_names[0]  # match argument to property name
                ret_names = ret_names[1:]
                lhs_len -= 1
            elif lhs_len > 0:
                lhs_kwargs[p.name] = lhs[lhs_idx]
        lhs_idx += 1

    # merge with existing
    lhs_kwargs.update(kwargs)
    return lhs_kwargs


def set_properties(algm, *args, **kwargs):
    """
        Set all of the properties of the algorithm. There is no guarantee of
        the order the properties will be set
        :param algm: An initialised algorithm object
        :param args: Positional arguments
        :param kwargs: Keyword arguments
    """
    def do_set_property(name, new_value):
        if new_value is None:
            return
        if isinstance(new_value, DataItem) and new_value.name():
            algm.setPropertyValue(name, new_value.name())
        else:
            algm.setProperty(name, new_value)
    # end
    if len(args) > 0:
        mandatory_props = algm.mandatoryProperties()
    else:
        mandatory_props = []

    postponed = []
    for (key, value) in iteritems(kwargs):
        if key in mandatory_props:
            mandatory_props.remove(key)
        if "IndexSet" in key:
            # The `IndexSet` sub-property of the "workspace property with index"
            # must be set after the workspace since it is validated based on in.
            postponed.append((key, value))
            continue
        do_set_property(key, value)
    for (key, value) in postponed:
        do_set_property(key, value)

    # zip stops at the length of the shorter list
    for (key, value) in zip(mandatory_props, args):
        do_set_property(key, value)


def is_function_property(prop):
    """
    Returns True if the property is a fit function

    :param prop: A property object
    :type Property
    :return:  True if the property is considered a fit function
    """
    return isinstance(prop, FunctionProperty)


def is_workspace_property(prop):

    """
        Returns true if the property is a workspace property.

        Currently several properties , i.e WorspaceProperty<EventWorkspace>
        cannot be recognised by Python so we have to resort to a name test

        :param prop: A property object
        :type Property
        :returns: True if the property is considered to be of type workspace
    """
    return isinstance(prop, IWorkspaceProperty)


def get_args_from_lhs(lhs, algm_obj):
    """
        Return the extra arguments that are to be passed to the algorithm
        from the information in the lhs tuple. These are basically the names
        of output workspaces.
        The algorithm properties are iterated over in the same order
        they were created within the wrapper and for each output
        workspace property an entry is added to the returned dictionary
        that contains {PropertyName:lhs_name}.

        :param lhs: A 2-tuple that contains the number of variables supplied on the lhs of the
        function call and the names of these variables
        :param algm_obj: An initialised algorithm object
        :returns: A dictionary mapping property names to the values extracted from the lhs variables
    """
    ret_names = lhs[1]
    extra_args = {}

    output_props = [algm_obj.getProperty(p) for p in algm_obj.outputProperties()]

    nprops = len(output_props)
    nnames = len(ret_names)

    name = 0

    for p in output_props:
        if is_workspace_property(p):
            # Check names is greater than 0 and less than nprops
            if 0 < nnames < nprops:
                extra_args[p.name] = ret_names[0]  # match argument to property name
                ret_names = ret_names[1:]
                nnames -= 1
            elif nnames > 0:
                extra_args[p.name] = ret_names[name]

        name += 1

    return extra_args


def _merge_keywords_with_lhs(keywords, lhs_args):
    """
        Merges the arguments from the two dictionaries specified
        by the keywords passed to a function and the lhs arguments
        that have been parsed. Any value in keywords overrides on
        in lhs_args.

        :param keywords: A dictionary of keywords that has been passed to the function call
        :param lhs_args: A dictionary of arguments retrieved from the lhs of the function call
    """
    final_keywords = lhs_args
    final_keywords.update(keywords)
    return final_keywords


def gather_returns(func_name, lhs, algm, ignore_regex=None, inout=False):
    """Gather the return values and ensure they are in the
       correct order as defined by the output properties and
       return them as a tuple. If their is a single return
       value it is returned on its own

       :param func_name: The name of the calling function.
       :param lhs: A 2-tuple that contains the number of variables supplied on the
       lhs of the function call and the names of these variables.
       :param algm: An executed algorithm object.
       :param ignore_regex: A list of strings containing regex expressions to match
       :param inout : gather also the InOut properties if True.
       against property names that will be ignored & not returned.
    """
    if ignore_regex is None:
        ignore_regex = []

    import re

    def ignore_property(name_to_check, regex_to_ignore):
        for regex in regex_to_ignore:
            if regex.match(name_to_check) is not None:
                return True
        # Matched nothing
        return False

    if type(ignore_regex) is str:
        ignore_regex = [ignore_regex]
    # Compile regexes
    for index, expr in enumerate(ignore_regex):
        ignore_regex[index] = re.compile(expr)

    retvals = OrderedDict()
    names = algm.outputProperties()
    if inout:
        names.extend(algm.inoutProperties())
    for name in names:
        if ignore_property(name, ignore_regex):
            continue
        prop = algm.getProperty(name)

        if is_workspace_property(prop):
            value = None
            if hasattr(prop, 'value'):
                value = prop.value
            if value is not None:
                retvals[name] = value
            else:
                try:
                    value_str = prop.valueAsStr
                    retvals[name] = AnalysisDataService.Instance()[value_str]
                except KeyError:
                    if not prop.isOptional() and prop.direction == Direction.InOut:
                        raise RuntimeError("Mandatory InOut workspace property '%s' on "
                                           "algorithm '%s' has not been set correctly. " % (name, algm.name()))
        elif is_function_property(prop):
            retvals[name] = FunctionWrapper(prop.value)
        else:
            if hasattr(prop, 'value'):
                retvals[name] = prop.value
            else:
                raise RuntimeError('Internal error. Unknown property type encountered. "%s" '
                                   'on algorithm "%s" is not understood by '
                                   'Python. Please contact development team' % (name, algm.name()))

    # If there is a snippet of code as follows
    # foo, bar, baz = simpleAPI.myFunc(...)
    # The number of values on LHS is 3 (foo, bar baz) and the number of
    # returned values is the number of values myFunc(...) returns
    number_of_returned_values = len(retvals)
    number_of_values_on_lhs = lhs[0]

    # If we have more than one value but not the same number of values throw
    if number_of_values_on_lhs > 1 and number_of_returned_values != number_of_values_on_lhs:
        # There is a discrepancy in the number are unpacking variables
        # Let's not have the more cryptic unpacking error raised
        raise RuntimeError("%s is trying to return %d output(s) but you have provided %d variable(s). "
                           "These numbers must match." % (func_name,
                                                          number_of_returned_values, number_of_values_on_lhs))
    if number_of_returned_values > 0:
        ret_type = namedtuple(func_name+"_returns", retvals.keys())
        ret_value = ret_type(**retvals)
        if number_of_returned_values == 1:
            return ret_value[0]
        else:
            return ret_value
    else:
        return None
