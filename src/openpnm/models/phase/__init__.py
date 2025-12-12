r"""
Phase
-----

This submodule contains models for calculating thermophysical properties of
liquids, gases, and solids.

"""

from ._phasedocs import *
from . import critical_props
from . import density
from . import diffusivity
from . import heat_capacity
from . import misc
from . import mixtures
from . import partition_coefficient
from . import surface_tension
from . import thermal_conductivity
from . import vapor_pressure
from . import viscosity


import logging
logger = logging.getLogger(__name__)


# %%
import inspect as _inspect
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
    'sigma': 'pore.surface_tension',
    'rhom': 'pore.molar_density',
    'rho': 'pore.density',
}


mixture_argmap = {
    'xs': 'pore.mole_fraction',
    'ys': 'pore.mole_fraction',
    'ks': 'pore.thermal_conductivity.*',
    'ws': 'pore.mass_fraction',
}


def chemicals_wrapper(phase, f, **kwargs):
    r"""
    Wrapper function for calling models in the ``chemicals`` package

    Parameters
    ----------
    phase : dict
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
        ``phase`` correspond to each argument required by ``f``. For
        instance, ``default_argmap['Tc'] = 'param.critical_temperature'``
        so any functions that need ``Tc`` will receive
        ``phase['param.critical_temperature']``. Some models only require
        the normal thermodynamic parameters, like ``T``, ``Tc``, and ``Pc``,
        so in many cases no additional arguments are needed beyond whats in the
        default argument map. However, other models require values that are
        themselves calcualted by another model, like ``Cvm``, or perhaps a
        list of constants like ``a0``, ``a1``, etc. It is necessary to provide
        these as keyword arguments, and they will be included in
        ``default_argmap`` as new entries or overwriting existing ones.
        So ``mu='pore.blah'`` will pass ``phase['pore.blah']`` to the
        ``mu`` argument of ``f``. Any arguments ending with an ``s`` are
        assumed to refer to the individual properties of pure components in a
        mixture. So ``mus`` means fetch the viscosity of each component as a
        list, like ``[1.3e-5, 1.9e-5]``. This function will trim the trailing
        ``s`` off any argument name before checking the ``default_argmap`` so
        that the normal pure component values can be looked up.

    Notes
    -----
    This wrapper works with both pure and mixture phases, but the mixture
    models are very slow due to the way ``chemicals`` vectorizes code. For pure
    species it allows the computation of values at many different conditions
    in a vectorized way. The means that the conditions in each pore, such as
    temperature, pressure, etc can be passed and iterpreted as a list of
    conditions. For mixture models, however, the vectorization is done over
    the compositions, at a *fixed* condition. This means that we must do a
    pure-python for-loop (ie. slow) for each individual pore. As such, we
    have re-implemented several of the most useful mixing models offered by
    ``chemicals`` in ``OpenPNM``, and include unit tests to ensure both
    implementations agree.

    """
    import chemicals as _chemicals
    # Update default argmap with any user supplied values
    argmap = default_argmap.copy()
    for k, v in kwargs.items():
        argmap[k] = v
    args = _get_items_from_target(phase, f, argmap)
    # f = getattr(_chemicals.numba_vectorized, f.__name__)
    msg = f"Numba version failed for {f.__name__}, reverting to pure python"
    if len(set(['xs', 'yz', 'zs']).intersection(args.keys())):
        # Call function in for-loop for each pore since they are not vectorized
        logger.info(msg)
        vals = _np.zeros(phase.Np)
        for pore in phase.Ps:
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
            f = getattr(_chemicals, f.__name__)
            vals = f(*args.values())
            logger.warn(msg)
    return vals


def _get_items_from_target(phase, f, argmap):
    import chemicals as _chemicals
    try:
        arginfo = _inspect.getfullargspec(f)
    except TypeError:
        temp = getattr(_chemicals, f.__name__)
        arginfo = _inspect.getfullargspec(temp)

    # Scan args and pull values from phase
    args = {}
    for item in arginfo.args:
        # Treat mole fraction specially, since it can be xs, ys, or zs
        if item in ['xs', 'ys', 'zs']:
            args[item] = list(phase['pore.mole_fraction'].values())
            continue
        if 'pore.' + item in phase.keys():
            # If eg. pore.Zc was added to dict directly, no argmap needed
            args[item] = phase[f"{'pore.'+item}"]
            continue
        if item in phase.params.keys():
            # If a parameter is in params but not in default_argmap
            args[item] = phase.params[item]
            continue
        if item.endswith('s'):  # Deal specifically with props that end in s
            if item[:-1] in argmap.keys():
                v = list(phase.get_comp_vals(argmap[item[:-1]]).values())
                args[item] = v
            elif item in argmap.keys():
                v = list(phase.get_comp_vals(argmap[item]).values())
                args[item] = v
        else:
            if item in argmap.keys():  # Get basic mixture properties like T&P
                args[item] = phase[argmap[item]]
    return args
