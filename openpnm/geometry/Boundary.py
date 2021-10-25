# -*- coding: utf-8 -*-
"""
===============================================================================
Boundary -- Subclass of GenericGeometry for Boundary Pores
===============================================================================

"""

from openpnm.models import geometry as gm
from openpnm.models import misc as mm
from openpnm.geometry import GenericGeometry


class Boundary(GenericGeometry):
    r"""
    Boundary subclass of GenericGeometry.

    The aim of this class is to assign geometric properties to boundary pores
    that are appropriate, such as 0 volume pores (for invasion simulations and
    porosity calculations) and very large throats (to offer negligible
    resistance in transport calculations).  When necessary, this class assumes
    that the pores are spherical and the throats are cylindrical.

    Parameters
    ----------
    network : OpenPNM Network object
        The Network to which the Geometry object should be associated

    pores : array_like
        The pores indice where the Geometry should be applied

    throats : array_like
        The throat indices where the Geometry should be applied

    name : string, optional
        A name to help identify the object, if not given one is automatically
        generated

    project : OpenPNM Project, optional
        Either a Network or a Project must be specified

    Examples
    --------
    >>> import openpnm
    >>> pn = openpnm.network.Cubic(shape=[3, 3, 3])
    >>> pn.add_boundary_pores()
    >>> Ps_int = pn.pores(labels=['*boundary'], mode='not')
    >>> Ps_boun = pn.pores(labels=['*boundary'])
    >>> Ts_int = pn.throats(labels=['*boundary'], mode='not')
    >>> Ts_boun = pn.throats(labels=['*boundary'])
    >>> geo = openpnm.geometry.StickAndBall(network=pn,
    ...                                     pores=Ps_int, throats=Ts_int)
    >>> boun = openpnm.geometry.Boundary(network=pn, pores=Ps_boun,
    ...                                  throats=Ts_boun)

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.add_model(propname='pore.area',
                       model=gm.pore_cross_sectional_area.sphere)

        self.add_model(propname='pore.volume',
                       model=mm.constant, value=0.0)

        self.add_model(propname='pore.seed',
                       model=mm.constant,
                       value=1.0)

        self.add_model(propname='throat.seed',
                       model=mm.constant, value=1.0)

        self.add_model(propname='pore.diameter',
                       model=mm.constant, value=0.0)

        self.add_model(propname='throat.diameter',
                       model=mm.from_neighbor_pores,
                       prop='pore.diameter', mode='max')

        self.add_model(propname='throat.endpoints',
                       model=gm.throat_endpoints.spherical_pores)

        self.add_model(propname='throat.length',
                       model=gm.throat_length.piecewise)

        self.add_model(propname='throat.conduit_lengths',
                       model=gm.throat_length.conduit_lengths)

        self.add_model(propname='throat.area',
                       model=gm.throat_cross_sectional_area.cylinder)

        self.add_model(propname='pore.area',
                       model=mm.from_neighbor_throats,
                       prop='throat.area', mode='max')

        self.add_model(propname='throat.surface_area',
                       model=gm.throat_surface_area.cylinder)

        self.add_model(propname='throat.volume',
                       model=gm.throat_volume.extrusion)
