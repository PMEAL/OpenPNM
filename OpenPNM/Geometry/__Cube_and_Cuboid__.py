# -*- coding: utf-8 -*-
"""
===============================================================================
Cube_and_Cuboid -- A standard Cubic pore and Cuboic throat model
===============================================================================

"""
import scipy as _sp
from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry import GenericGeometry


class Cube_and_Cuboid(GenericGeometry):
    r"""
    Toray090 subclass of GenericGeometry

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.models.add(propname='pore.seed',
                        model=gm.pore_misc.random)
        # Find Network spacing
        Ps = self._net.pores(self.name)
        Ts = self._net.find_neighbor_throats(pores=Ps, mode='intersection')
        P1 = self._net['throat.conns'][:, 0][Ts]
        P2 = self._net['throat.conns'][:, 1][Ts]
        C1 = self._net['pore.coords'][P1]
        C2 = self._net['pore.coords'][P2]
        E = _sp.sqrt(_sp.sum((C1-C2)**2, axis=1))  # Euclidean distance
        if _sp.allclose(E, E[0]):
            spacing = E[0]
        else:
            raise Exception('A unique value of spacing could not be inferred')
        self.models.add(propname='pore.diameter',
                        model=gm.pore_diameter.normal,
                        loc=spacing/2,
                        scale=spacing/10)
        self.models.add(propname='pore.area',
                        model=gm.pore_area.cubic)
        self.models.add(propname='pore.volume',
                        model=gm.pore_volume.cube)
        self.models.add(propname='throat.diameter',
                        model=gm.throat_diameter.minpore,
                        factor=0.5)
        self.models.add(propname='throat.length',
                        model=gm.throat_length.straight)
        self.models.add(propname='throat.volume',
                        model=gm.throat_volume.cuboid)
        self.models.add(propname='throat.area',
                        model=gm.throat_area.cuboid)
        self.models.add(propname='throat.surface_area',
                        model=gm.throat_surface_area.cuboid)
