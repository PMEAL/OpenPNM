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


from . import critical_props
from . import density
from . import diffusivity
from . import electrical_conductivity
from . import heat_capacity
from . import misc
from . import mixtures
from . import molar_density
from . import partition_coefficient
from . import permittivity
from . import surface_tension
from . import thermal_conductivity
from . import vapor_pressure
from . import viscosity


import inspect as _inspect
import chemicals as _chemicals
import numpy as _np


_argmap = {
    'T': 'pore.temperature',
    'P': 'pore.pressure',
    'MW': 'param.molecular_weight',
    'Tc': 'param.critical_temperature',
    'Pc': 'param.critical_pressure',
    'Zc': 'param.critical_compressibilty_factor',
    'Vc': 'param.critical_volume',
    'Tm': 'param.melting_temperature',
    'Tb': 'param.boiling_temperature',
    'omega': 'param.acentric_factor',
    'dipole_moment': 'param.dipole_moment',
}


def chemicals_pure_prop(target, f, argmap=None):
    if argmap is None:
        argmap = _argmap
    # Get the non-numba version of f, which is needed for inspect to extract info
    temp = getattr(_chemicals, f.__name__)
    arginfo = _inspect.getfullargspec(temp)
    # Scan args and set defaults for any that are NOT in argmap
    if arginfo.defaults is not None:
        offset = len(arginfo.args) - len(arginfo.defaults)
        for i in _np.arange(len(arginfo.defaults)):
            if arginfo.args[offset + i] not in argmap.keys():
                argmap[arginfo.args[offset + i]] = arginfo.defaults[i]
    # Scan through the arguments and get necessary values from target
    args = {}
    for item in arginfo.args:
        if argmap[item] is None:
            args[item] = None
        elif argmap[item] == '':
            args[item] = ''
        else:
            args[item] = target[argmap[item]]
    try:
        # Get the numba vectorized version of f, or else numpy arrays don't work
        f = getattr(_chemicals.numba_vectorized, f.__name__)
        # Compute values
        vals = f(*args.values())
    except:
        raise Exception('numba vectorized version did not work')
    return vals
