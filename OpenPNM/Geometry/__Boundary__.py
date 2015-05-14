# -*- coding: utf-8 -*-
"""
===============================================================================
Boundary -- Subclass of GenericGeometry for Boundary Pores
===============================================================================

"""

from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry import GenericGeometry


class Boundary(GenericGeometry):
    r"""
    Boundary subclass of GenericGeometry.

    Parameters
    ----------
    network : OpenPNM Network object
        The Network to which the Geometry object should be associated
    pores, throats : array_like
        The pores and/or throats where the Geometry should be applied
    shape: str
        Stick and Ball or Cube and Cuboid? ('spheres','cubes')

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> Ps_int = pn.pores(labels=['top', 'bottom'], mode='not')
    >>> Ps_boun = pn.pores(labels=['top', 'bottom'], mode='union')
    >>> Ts_int = pn.throats(labels=['top', 'bottom'], mode='not')
    >>> Ts_boun = pn.throats(labels=['top', 'bottom'], mode='union')
    >>> geo = OpenPNM.Geometry.Cube_and_Cuboid(network=pn,
    ...                                        pores=Ps_int,
    ...                                        throats=Ts_int)
    >>> boun = OpenPNM.Geometry.Boundary(network=pn, pores=Ps_boun, throats=Ts_boun)

    """

    def __init__(self, shape='spheres', **kwargs):
        super().__init__(**kwargs)
        self._generate(shape)

    def _generate(self, shape):
        try:
            self['pore.seed']
            seeds = True
        except:
            seeds = False

        if seeds:
            self.models.add(propname='pore.seed',
                            model=gm.pore_misc.constant,
                            value=0.9999)
        self.models.add(propname='pore.diameter',
                        model=gm.pore_misc.constant,
                        value=0)
        if seeds:
            self.models.add(propname='throat.seed',
                            model=gm.throat_misc.neighbor,
                            pore_prop='pore.seed',
                            mode='max')
        self.models.add(propname='throat.diameter',
                        model=gm.throat_misc.neighbor,
                        pore_prop='pore.diameter',
                        mode='max')
        self['pore.volume'] = 0.0
        self['pore.seed'] = 0.0
        self.models.add(propname='throat.length', model=gm.throat_length.straight)
        self['throat.volume'] = 0.0
        self['throat.seed'] = 0.0
        if shape == 'spheres':
            self.models.add(propname='throat.area', model=gm.throat_area.cylinder)
            self.models.add(propname='throat.surface_area',
                            model=gm.throat_surface_area.cylinder)
        elif shape == 'cubes':
            self.models.add(propname='throat.area', model=gm.throat_area.cuboid)
            self.models.add(propname='throat.surface_area',
                            model=gm.throat_surface_area.cuboid)
        self['pore.area'] = 1.0
