# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 09:26:18 2018

@author: Precision
"""

import openpnm as op
import numpy as np
import matplotlib.pyplot as plt
import openpnm.utils.vertexops as vo
from openpnm.geometry import models as gm
plt.close('all')

np.random.seed(1)
scale = 5e-5
coords = np.random.random([50, 3])*scale
ws = op.core.Workspace()
ws.loglevel = 20
sim = op.core.Simulation()
sh = [scale, scale, scale]
dual_net = op.materials.VoronoiFibers(points=coords,
                                      shape=sh,
                                      name='dual',
                                      simulation=sim,
                                      fiber_rad=2e-6,
                                      vox_len=1e-6)
dual_del = sim.geometries['dual_del']
dual_del.plot_porosity_profile()
plt.figure()
plt.imshow(dual_del._fiber_image[:, :, 50])
plt.figure()
plt.imshow(dual_del._hull_image[:, :, 50])
plt.figure()
plt.imshow(dual_del._dt_image[:, :, 50])

vo.plot_throat(dual_del, [0])

dual_vor = sim.geometries['dual_vor']

