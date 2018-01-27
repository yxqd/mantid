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


@fitting_algorithm(inout=True)
def Fit(*args, **kwargs):
    """
    Fit defines the interface to the fitting within Mantid.
    It can work with arbitrary data sources and therefore some options
    are only available when the function & workspace type are known.

    This simple wrapper takes the Function (as a string or a
    FunctionWrapper object) and the InputWorkspace
    as the first two arguments. The remaining arguments must be
    specified by keyword.

    Example:
      Fit(Function='name=LinearBackground,A0=0.3', InputWorkspace=dataWS',
          StartX='0.05',EndX='1.0',Output="Z1")
    """
    return None
