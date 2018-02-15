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
        self.models = self.recipe()
        self.regenerate_models()

    @classmethod
    def recipe(cls):
        r = {'pore.seed': {'model': gm.pore_misc.random,
                           'num_range': [0, 0.1],
                           'seed': None},

             'pore.max_size': {'model': gm.pore_size.largest_sphere,
                               'iters': 10,
                               'pore_diameter': 'pore.diameter'},

             'pore.volume': {'model': gm.pore_volume.sphere,
                             'pore_diameter': 'pore.diameter'},

             'pore.diameter': {'model': gm.misc.product,
                               'arg1': 'pore.max_size',
                               'arg2': 'pore.seed'},
                               
             'pore.area': {'model': gm.pore_area.spherical,
                             'pore_diameter': 'pore.diameter'},                               
                               
             'throat.max_size': {'model': gm.throat_misc.neighbor,
                                 'mode': 'min',
                                 'pore_prop': 'pore.diameter'},

             'throat.diameter': {'model': gm.misc.scaled,
                                 'factor': 0.5,
                                 'prop': 'throat.max_size'},

             'throat.length': {'model': gm.throat_length.straight,
                               'L_negative': 1e-09,
                               'pore_diameter': 'pore.diameter'},

             'throat.surface_area': {'model': gm.throat_surface_area.cylinder,
                                     'throat_diameter': 'throat.diameter',
                                     'throat_length': 'throat.length'},

             'throat.volume': {'model': gm.throat_volume.cylinder,
                               'throat_diameter': 'throat.diameter',
                               'throat_length': 'throat.length'},

             'throat.area': {'model': gm.throat_area.cylinder,
                             'throat_diameter': 'throat.diameter'},

             }
        return r
