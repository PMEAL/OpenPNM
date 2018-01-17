# -*- coding: utf-8 -*-
"""
===============================================================================
Stick_and_Ball -- A standard 'stick & ball' geometrical model
===============================================================================

"""
from openpnm.core import ModelsDict
from openpnm.geometry import models as gm
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

    def __init__(self, network, **kwargs):
        super().__init__(network, **kwargs)

        self.add_model(propname='pore.max_size',
                       model=gm.pore_size.largest_sphere)
        self.add_model(propname='pore.seed',
                       model=gm.pore_misc.random,
                       seed=None, num_range=[0, 0.1],
                       regen_mode='constant')
        self.add_model(propname='pore.diameter',
                       model=gm.misc.product,
                       arg1='pore.max_size', arg2='pore.seed')
        self.add_model(propname='pore.volume',
                       model=gm.pore_volume.sphere)
        self.add_model(propname='throat.max_size',
                       model=gm.throat_misc.neighbor,
                       pore_prop='pore.diameter', mode='min')
        self.add_model(propname='throat.diameter',
                       model=gm.misc.scaled,
                       prop='throat.max_size', factor=0.5)
        self.add_model(propname='throat.length',
                       model=gm.throat_length.straight)
        self.add_model(propname='throat.volume',
                       model=gm.throat_volume.cylinder)
        self.add_model(propname='throat.area',
                       model=gm.throat_area.cylinder)
        self.add_model(propname='throat.surface_area',
                       model=gm.throat_surface_area.cylinder)
#        self.models.update(self.recipe())
        self.regenerate_models()

    @classmethod
    def recipe(cls):
        r = {'pore.seed': {'model': gm.pore_misc.random,
                           'num_range': [0, 0.1],
                           'seed': None},

             'pore.max_size': {'model': gm.pore_size.largest_sphere,
                               'iters': 10,
                               'pore_diameter': 'pore.diameter'},

             'pore.diameter': {'model': gm.misc.product,
                               'arg1': 'pore.max_size',
                               'arg2': 'pore.seed'},

             'pore.volume': {'model': gm.pore_volume.sphere,
                             'pore_diameter': 'pore.diameter'},

             'throat.area': {'model': gm.throat_area.cylinder,
                             'throat_diameter': 'throat.diameter'},

             'throat.diameter': {'model': gm.misc.scaled,
                                 'factor': 0.5,
                                 'prop': 'throat.max_size'},

             'throat.length': {'model': gm.throat_length.straight,
                               'L_negative': 1e-09,
                               'pore_diameter': 'pore.diameter'},

             'throat.max_size': {'model': gm.throat_misc.neighbor,
                                 'mode': 'min',
                                 'pore_prop': 'pore.diameter'},

             'throat.surface_area': {'model': gm.throat_surface_area.cylinder,
                                     'throat_diameter': 'throat.diameter',
                                     'throat_length': 'throat.length'},

             'throat.volume': {'model': gm.throat_volume.cylinder,
                               'throat_diameter': 'throat.diameter',
                               'throat_length': 'throat.length'}}
        return r
