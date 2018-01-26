# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 09:26:18 2018

@author: Precision
"""

import openpnm as op
import numpy as np
import matplotlib.pyplot as plt
import openpnm.utils.vertexops as vo
plt.close('all')

np.random.seed(1)
scale = 1e-4
coords = np.random.random([100, 3])*scale
ws = op.core.Workspace()
ws.loglevel = 20
sim = op.core.Simulation()
sh = [scale, scale, scale]
vor = op.materials.VoronoiFibers(points=coords, shape=sh, name='dual',
                                 simulation=sim,
                                 fiber_rad=3e-6,
                                 vox_len=1e-6)
geom = sim.geometries['dual_del']
geom.plot_porosity_profile()
plt.figure()
plt.imshow(geom._fiber_image[:, :, 50])
plt.figure()
plt.imshow(geom._hull_image[:, :, 50])
plt.figure()
plt.imshow(geom._dt_image[:, :, 50])

vo.plot_throat(geom, [0])
