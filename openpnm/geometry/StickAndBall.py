# -*- coding: utf-8 -*-
"""
===============================================================================
Stick_and_Ball -- A standard 'stick & ball' geometrical model
===============================================================================

"""
import scipy as _sp
from openpnm.core import ModelsMixin
from openpnm.geometry import models as gm
from openpnm.geometry import GenericGeometry


class StickAndBall(GenericGeometry):
    r"""
    Stick and Ball subclass of GenericGeometry.  This subclass is meant as a
    basic default geometry to get started quickly.

    Parameters
    ----------
    name : string
        The name of the object, which is also used as the label where this
        geometry is defined.

    """

    def __init__(self, network, **kwargs):
        super().__init__(network, **kwargs)

        # Find Network spacing
        Ps = network.pores(self.name)
        Ts = network.find_neighbor_throats(pores=Ps, mode='intersection')
        P1 = network['throat.conns'][:, 0][Ts]
        P2 = network['throat.conns'][:, 1][Ts]
        C1 = network['pore.coords'][P1]
        C2 = network['pore.coords'][P2]
        E = _sp.sqrt(_sp.sum((C1-C2)**2, axis=1))  # Euclidean distance
        if _sp.allclose(E, E[0]):
            spacing = E[0]
        else:
            raise Exception('A unique value of spacing could not be inferred')

        self.add_model(propname='pore.seed',
                       model=gm.pore_misc.random,
                       seed=1, num_range=[0, 0.1],
                       regen_mode='constant')
        self.add_model(propname='pore.diameter',
                       model=gm.pore_size.normal,
                       loc=spacing/2,
                       scale=spacing/10)
        self.add_model(propname='pore.volume',
                       model=gm.pore_volume.sphere)
        self.add_model(propname='throat.diameter',
                       model=gm.throat_misc.neighbor,
                       pore_prop='pore.diameter', mode='min')
        self.add_model(propname='throat.length',
                       model=gm.throat_length.straight)
        self.add_model(propname='throat.volume',
                       model=gm.throat_volume.cylinder)
        self.add_model(propname='throat.area',
                       model=gm.throat_area.cylinder)
        self.add_model(propname='throat.surface_area',
                       model=gm.throat_surface_area.cylinder)
        self.regenerate_models()
