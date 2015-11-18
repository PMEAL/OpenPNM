# -*- coding: utf-8 -*-
"""
===============================================================================
PeriodicThroats -- A basic Geometry
===============================================================================

"""

import scipy as _sp
from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry import GenericGeometry


class PeriodicThroats(GenericGeometry):
    r"""
    This class is meant to calculate geometric properties for periodic throats.
    The properties will be such that when regular pore scale models are applied
    they will be fed with correct values to produce the desired results. The
    main issue with these throats is the length, which if calculated using the
    normal way would be as long as the domain.

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network to which this Geometry object should be attached

    throats : array_like
        The throats to which this Geometry object should be applied

    name : string
        The name to assign to the object.  This can be helpful identifying
        objects.  If no name is sent a random one will be generated.

    Notes
    -----
    It is not possible to assign this PeriodicThroat Geometry to any pores, so
    an exception will be rasied if any pores are recieved as arguments.

    """

    def __init__(self, **kwargs):
        if 'pores' in kwargs:
            raise Exception('This Geometry subclass is intended for throats only')

        super().__init__(**kwargs)

        self['throat.length'] = 0
        self['throat.diameter'] = _sp.inf
