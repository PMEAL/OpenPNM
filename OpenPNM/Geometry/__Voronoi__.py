# -*- coding: utf-8 -*-
"""
===============================================================================
Voronoi --Subclass of GenericGeometry for a standard Geometry created from a
Voronoi Diagram Used with Delaunay Network but could work for others (not tested)
===============================================================================


"""

import scipy as sp
from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry import GenericGeometry


class Voronoi(GenericGeometry):
    r"""
    Voronoi subclass of GenericGeometry.

    Parameters
    ----------

    """

    def __init__(self, fibre_rad=3e-06, load_gen='gen', **kwargs):
        super().__init__(**kwargs)
        if load_gen == 'gen':
            self._generate(fibre_rad)

    def _generate(self, fibre_rad):
        self.models.add(propname='pore.vertices',
                        model=gm.pore_vertices.voronoi)
        self.models.add(propname='throat.vertices',
                        model=gm.throat_vertices.voronoi)
        self.models.add(propname='throat.normal',
                        model=gm.throat_normal.voronoi)
        self.models.add(propname='throat.offset_vertices',
                        model=gm.throat_offset_vertices.distance_transform,
                        offset=fibre_rad,
                        set_dependent=True)
        # Remove throats that are fully occluded and pores with no connections
        self._net.trim_occluded_throats()

        self.models.add(propname='pore.seed',
                        model=gm.pore_misc.random,
                        seed=self._seed)
        self.models.add(propname='throat.seed',
                        model=gm.throat_misc.neighbor,
                        pore_prop='pore.seed',
                        mode='min')
        self.models.add(propname='pore.volume',
                        model=gm.pore_volume.voronoi)
        self.models.add(propname='pore.diameter',
                        model=gm.pore_diameter.voronoi)
        self.models.add(propname='pore.area',
                        model=gm.pore_area.spherical)
        self.models.add(propname='throat.diameter',
                        model=gm.throat_diameter.voronoi)
        self.models.add(propname='throat.length',
                        model=gm.throat_length.constant,
                        const=fibre_rad*2)
        self.models.add(propname='throat.volume',
                        model=gm.throat_volume.extrusion)
        self.models.add(propname='throat.surface_area',
                        model=gm.throat_surface_area.extrusion)
