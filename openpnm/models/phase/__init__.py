r"""
Phase
-----

This submodule contains models for calculating thermophysical properties of
liquids, gases, and solids.

"""

# The following bits are to initialize some boilerplate docstrings for docrep
from openpnm.utils import Docorator as _doc
_docstr = _doc()
_docstr.params['models.phase.T'] = \
    r"""temperature : str
            Name of the dictionary key on ``target`` where the array containing
            temperature values is stored"""
_docstr.params['models.phase.P'] = \
    r"""pressure : str
            Name of the dictionary key on ``target`` where the array containing
            pressure values is stored"""
_docstr.params['models.phase.SI_note'] = \
    r"""Since numpy arrays don't natively support units, OpenPNM cannot enforce
    units on values and arrays, however the use of SI is assumed."""


from . import density
from . import diffusivity
from . import electrical_conductivity
from . import misc
from . import mixtures
from . import molar_density
from . import partition_coefficient
from . import permittivity
from . import surface_tension
from . import thermal_conductivity
from . import vapor_pressure
from . import viscosity
