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

# local imports
from .fit_common import fitting_algorithm


@fitting_algorithm()
def EvaluateFunction(*args, **kwargs):
    """
    This function evaluates a function on a data set.
    The data set is defined in a way similar to Fit algorithm.

    Example:
      EvaluateFunction(Function='name=LinearBackground,A0=0.3', InputWorkspace=dataWS',
          StartX='0.05',EndX='1.0',Output="Z1")
    """
    return None
