# -*- coding: utf-8 -*-
import openpnm.models.geometry as gm
from openpnm.geometry import GenericGeometry


class StickAndBall(GenericGeometry):
    r"""
    Stick and Ball subclass of GenericGeometry.  This subclass is meant as a
    basic default geometry to get started quickly.

    Pore diameters are randomly assigned between 0 and the largest sphere that
    does not overlap with it's nearest neighbor.

    Throat diameters are half the diameter of the smaller of it's two
    neighboring pores.

    Parameters
    ----------
    name : string
        The name of the object, which is also used as the label where this
        geometry is defined.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.add_model(propname='pore.seed',
                       model=gm.pore_misc.random,
                       num_range=[0, 0.1],
                       seed=None)

        self.add_model(propname='pore.max_size',
                       model=gm.pore_size.largest_sphere,
                       iters=10)

        self.add_model(propname='pore.area',
                       model=gm.pore_area.spherical,
                       pore_diameter='pore.diameter')

        self.add_model(propname='pore.volume',
                       model=gm.pore_volume.sphere,
                       pore_diameter='pore.diameter')

        self.add_model(propname='pore.diameter',
                       model=gm.misc.product,
                       arg1='pore.max_size',
                       arg2='pore.seed')

        self.add_model(propname='throat.max_size',
                       model=gm.throat_misc.neighbor,
                       mode='min',
                       pore_prop='pore.diameter')

        self.add_model(propname='throat.diameter',
                       model=gm.misc.scaled,
                       factor=0.5,
                       prop='throat.max_size')

        self.add_model(propname='throat.length',
                       model=gm.throat_length.straight,
                       L_negative=1e-09,
                       pore_diameter='pore.diameter')

        self.add_model(propname='throat.surface_area',
                       model=gm.throat_surface_area.cylinder,
                       throat_diameter='throat.diameter',
                       throat_length='throat.length')

        self.add_model(propname='throat.volume',
                       model=gm.throat_volume.cylinder,
                       throat_diameter='throat.diameter',
                       throat_length='throat.length')

        self.add_model(propname='throat.area',
                       model=gm.throat_area.cylinder,
                       throat_diameter='throat.diameter')
