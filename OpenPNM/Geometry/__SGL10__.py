# -*- coding: utf-8 -*-
"""
===============================================================================
SGL10 -- A geometry model for SGL10 type Gas Diffusion Layers
===============================================================================

"""

from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry import GenericGeometry


class SGL10(GenericGeometry):
    r"""
    SGL10 subclass of GenericGeometry.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._generate()

    def _generate(self):
        self.models.add(propname='pore.seed',
                        model=gm.pore_misc.random,
                        num_range=[0, 0.8834],
                        seed=self._seed,
                        regen_mode='constant')
        self.models.add(propname='throat.seed',
                        model=gm.throat_misc.neighbor,
                        pore_prop='pore.seed',
                        mode='min')
        self.models.add(propname='pore.diameter',
                        model=gm.pore_diameter.sphere,
                        psd_name='weibull_min',
                        psd_shape=3.07,
                        psd_loc=1.97e-6,
                        psd_scale=1.6e-5,
                        psd_offset=18e-6)
        self.models.add(propname='pore.area',
                        model=gm.pore_area.spherical)
        self.models.add(propname='pore.volume',
                        model=gm.pore_volume.sphere)
        self.models.add(propname='throat.diameter',
                        model=gm.throat_diameter.cylinder,
                        tsd_name='weibull_min',
                        tsd_shape=3.07,
                        tsd_loc=1.97e-6,
                        tsd_scale=1.6e-5,
                        tsd_offset=18e-6)
        self.models.add(propname='throat.length',
                        model=gm.throat_length.straight)
        self.models.add(propname='throat.volume',
                        model=gm.throat_volume.cylinder)
        self.models.add(propname='throat.area',
                        model=gm.throat_area.cylinder)
        self.models.add(propname='throat.surface_area',
                        model=gm.throat_surface_area.cylinder)
