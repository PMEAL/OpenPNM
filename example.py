# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 22:51:50 2019

@author: Jeff
"""

import openpnm as op
print(op.__version__)
from openpnm.phases import mixtures

pn = op.network.Cubic(shape=[3, 3])
geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)

N2 = mixtures.species.gases.N2(network=pn, name='nitrogen')
O2 = mixtures.species.gases.O2(network=pn, name='oxygen')

air = mixtures.GenericMixture(network=pn, components=[N2, O2])
air.set_mole_fraction(component=O2, values=0.21)
air.set_mole_fraction(component=N2, values=0.79)
air.add_model(propname='pore.diffusivity',
              model=op.models.phases.mixtures.wilke_fuller_diffusivity)

phys = op.physics.Standard(network=pn, phase=air, geometry=geo)


