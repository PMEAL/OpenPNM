# -*- coding: utf-8 -*-
"""
===============================================================================
Stick_and_Ball -- A standard 'stick & ball' geometrical model
===============================================================================

"""

from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry import GenericGeometry


class Stick_and_Ball(GenericGeometry):
    r"""
    Stick and Ball subclass of GenericGeometry.  This subclass is meant as a
    basic default geometry to get started quickly.

    Parameters
    ----------
    name : string
        The name of the object, which is also used as the label where this
        geometry is defined.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._generate()

    def _generate(self):
        self.models.add(propname='pore.seed',
                        model=gm.pore_misc.random,
                        regen_mode='constant')
        self.models.add(propname='pore.diameter',
                        model=gm.pore_diameter.normal,
                        loc=self._net._spacing[0]/2,
                        scale=self._net._spacing[0]/10)
        self.models.add(propname='pore.area',
                        model=gm.pore_area.spherical)
        self.models.add(propname='pore.volume',
                        model=gm.pore_volume.sphere)
        self.models.add(propname='throat.diameter',
                        model=gm.throat_diameter.minpore,
                        factor=0.5)
        self.models.add(propname='throat.length',
                        model=gm.throat_length.straight)
        self.models.add(propname='throat.volume',
                        model=gm.throat_volume.cylinder)
        self.models.add(propname='throat.area',
                        model=gm.throat_area.cylinder)
        self.models.add(propname='throat.surface_area',
                        model=gm.throat_surface_area.cylinder)
