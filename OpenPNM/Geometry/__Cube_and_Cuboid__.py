# -*- coding: utf-8 -*-
"""
===============================================================================
Cube_and_Cuboid -- A standard Cubic pore and Cuboic throat model
===============================================================================

"""

from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry import GenericGeometry


class Cube_and_Cuboid(GenericGeometry):
    r"""
    Toray090 subclass of GenericGeometry

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._generate()

    def _generate(self):
        self.models.add(propname='pore.seed',
                        model=gm.pore_misc.random)
        self.models.add(propname='pore.diameter',
                        model=gm.pore_diameter.normal,
                        loc=self._net._spacing[0]/2,
                        scale=self._net._spacing[0]/10)
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
