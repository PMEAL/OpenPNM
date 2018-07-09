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

    Parameters
    ----------
    network : OpenPNM Network object
        The Network to which the Geometry object should be associated

    pores : array_like
        The pores indice where the Geometry should be applied

    throat : array_like
        The throat indices where the Geometry should be applied

    shape: str
        Stick and Ball or Cube and Cuboid? ('spheres','cubes')

    name : string, optional
        A name to help identify the object, if not given one is automatically
        generated

    project : OpenPNM Project, optional
        Either a Network or a Project must be specified

    Examples
    --------
    >>> import openpnm
    >>> pn = openpnm.network.Cubic(shape=[3, 3, 3])
    >>> Ps_int = pn.pores(labels=['top', 'bottom'], mode='not')
    >>> Ps_boun = pn.pores(labels=['top', 'bottom'], mode='union')
    >>> Ts_int = pn.throats(labels=['internal'])
    >>> Ts_boun = pn.throats(labels=['internal'], mode='not')
    >>> geo = openpnm.geometry.GenericGeometry(network=pn,
    ...                                        pores=Ps_int, throats=Ts_int)
    >>> boun = openpnm.geometry.Boundary(network=pn, pores=Ps_boun,
    ...                                  throats=Ts_boun)

    """

    def __init__(self, shape='spheres', **kwargs):
        super().__init__(**kwargs)
        self['pore.diameter'] = 1e-12
        self.add_model(propname='throat.diameter',
                       model=mm.from_neighbor_pores,
                       pore_prop='pore.diameter',
                       mode='max')
        self['pore.volume'] = 1e-36
        self['pore.seed'] = 1.0
        self['throat.seed'] = 1.0
        self['throat.volume'] = 1e-36
        self.add_model(propname='throat.length',
                       model=gm.throat_length.straight)
        if shape == 'spheres':
            self.add_model(propname='throat.area',
                           model=gm.throat_area.cylinder)
            self.add_model(propname='throat.surface_area',
                           model=gm.throat_surface_area.cylinder)
        elif shape == 'cubes':
            self.add_model(propname='throat.area',
                           model=gm.throat_area.cuboid)
            self.add_model(propname='throat.surface_area',
                            model=gm.throat_surface_area.cuboid)
