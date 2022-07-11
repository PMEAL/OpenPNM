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
    'dipole': 'param.dipole_moment',
    'mu': 'pore.viscosity',
}


def chemicals_wrapper(target, f, **kwargs):
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
        ``mu`` argument of ``f``. Any arguments ending with an ``s`` are
        assumed to refer to the individual properties of pure components in a
        mixture. So ``mus`` means fetch the viscosity of each component as a
        list, like ``[1.3e-5, 1.9e-5]``. This function will trim the trailing
        ``s`` off any argument name before checking the ``default_argmap`` so
        that the normal pure component values can be looked up.

    """
    # Update default argmap with any user supplied values
    argmap = default_argmap.copy()
    for k, v in kwargs.items():
        argmap[k] = v
    args = _get_items_from_target(target, f, argmap)
    # f = getattr(_chemicals.numba_vectorized, f.__name__)
    if any([True for k in args.keys() if k.endswith('s')]):
        # Call function in for-loop for each pores since they are not vectorized
        vals = _np.zeros(target.Np)
        for pore in target.Ps:
            a = {}
            for item in args.keys():
                if item.endswith('s'):
                    try:
                        a[item] = [args[item][i][pore] for i in range(len(args[item]))]
                    except TypeError:
                        a[item] = args[item]
                else:
                    a[item] = args[item][pore]
            vals[pore] = f(*list(a.values()))
    else:
        try:
            # Get the numba vectorized version of f, or else numpy arrays don't work
            f = getattr(_chemicals.numba_vectorized, f.__name__)
            vals = f(*args.values())
        except AssertionError:
            raise Exception('numba vectorized version did not work')
    return vals


def _get_items_from_target(target, f, argmap):
    try:
        arginfo = _inspect.getfullargspec(f)
    except TypeError:
        temp = getattr(_chemicals, f.__name__)
        arginfo = _inspect.getfullargspec(temp)

    # Scan args and pull values from target
    args = {}
    for item in arginfo.args:
        # Treat mole fraction specially, since it can be xs, ys, or zs
        if item in ['xs', 'ys', 'zs']:
            args[item] = list(target['pore.mole_fraction'].values())
        elif item in argmap.keys():  # Get basic mixture properties like T&P
            args[item] = target[argmap[item]]
        elif item[:-1] in argmap.keys():  # Get props that end in s
            v = list(target.get_comp_vals(argmap[item[:-1]]).values())
            args[item] = v
        elif 'pore.' + item in target.keys():
            # If eg. pore.Zc was added to dict directly, no argmap needed
            args[item] = target[f"{'pore.'+item}"]
        elif item in target.params.keys():
            # If a parameter is in params but not in default_argmap
            args[item] = target.params[item]
    return args
