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


default_argmap = {
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
    'mu': 'pore.viscosity',
}


def chemicals_pure_prop(target, f, **kwargs):
    r"""
    Wrapper function for calling models in the ``chemicals`` package

    Parameters
    ----------
    target : dict
        The OpenPNM Species object for which this model should calculate
        values. This object should ideally have all the necessary chemical
        properties in its ``params`` attribute, although this is optional as
        discussed below in relation to ``kwargs``.
    f : function
        The handle of the function to apply, such as
        ``chemicals.viscosity.Letsou_Stiel``.
    kwargs
        By default this function will use
        ``openpnm.models.phase.default_argmap`` to determine which names on
        ``target`` correspond to each argument required by ``f``. For
        instance, ``default_argmap['Tc'] = 'param.critical_temperature'``
        so any functions that need ``Tc`` will receive
        ``target['param.critical_temperature']``. Some models only require
        the normal thermodynamic parameters, like ``T``, ``Tc``, and ``Pc``,
        so in many cases no additional arguments are needed. However, other
        models require values that are themselves calcualted by another model,
        like ``Cvm``, or perhaps a list of constants like ``a0``, ``a1``, etc.
        It is necessary to provide these as keyword arguments, and they will be
        included in ``default_argmap`` as new entries or overwriting existing
        ones. So ``mu='pore.blah'`` will pass ``target['pore.blah']`` to the
        ``mu`` argument of ``f``.

    """
    # Update default argmap with any user supplied values
    for k, v in kwargs.items():
        default_argmap['k'] = v
    # Get the non-numba version of f, which is needed for inspect to extract info
    temp = getattr(_chemicals, f.__name__)
    arginfo = _inspect.getfullargspec(temp)
    # Scan args and set defaults for any that are NOT in argmap
    if arginfo.defaults is not None:
        offset = len(arginfo.args) - len(arginfo.defaults)
        for i in _np.arange(len(arginfo.defaults)):
            if arginfo.args[offset + i] not in default_argmap.keys():
                default_argmap[arginfo.args[offset + i]] = arginfo.defaults[i]
    # Scan through the arguments and get necessary values from target
    args = {}
    for item in arginfo.args:
        # Map the parameter (eg Tc) to dict name (eg param.critical_temperature)
        if item in default_argmap.keys():
            if default_argmap[item] is None:
                args[item] = None
            elif default_argmap[item] == '':
                args[item] = ''
            else:
                args[item] = target[default_argmap[item]]
        else:
            # If item is not in argmap, then just try getting it on target
            args[item] = target['pore.'+item]
    try:
        # Get the numba vectorized version of f, or else numpy arrays don't work
        f = getattr(_chemicals.numba_vectorized, f.__name__)
        # Compute values
        vals = f(*args.values())
    except:
        raise Exception('numba vectorized version did not work')
    return vals
