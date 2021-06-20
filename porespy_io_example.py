import openpnm as op
import porespy as ps
import matplotlib.pyplot as plt
import numpy as np

im = ps.generators.blobs([300, 300, 300], porosity=0.75)
im = ps.filters.fill_blind_pores(im)
# plt.imshow(im)

snow = ps.networks.snow2(im)
# plt.imshow(snow.regions/snow.phases)

# %%
pn, geo = op.io.PoreSpy.import_data(snow.network,
                                    settings={'pore_shape': 'cone',
                                              'throat_shape': 'cylinder'})
h = pn.check_network_health()
op.topotools.trim(network=pn, pores=h['trim_pores'])
air = op.phases.Air(network=pn)
phys = op.physics.Basic(network=pn, phase=air, geometry=geo)

fd = op.algorithms.FickianDiffusion(network=pn, phase=air)
fd.set_value_BC(pores=pn.pores('ymin'), values=1.0)
fd.set_value_BC(pores=pn.pores('ymax'), values=0.0)
fd.run()

# %%  Plot the concentration map over the image
conc = fd['pore.concentration'][snow.regions - 1]*snow.phases
fig = op.topotools.plot_connections(pn)
fig = op.topotools.plot_coordinates(pn, fig=fig)
plt.imshow((conc/snow.phases).T, origin='upper', interpolation='none')
plt.axis(False)
Deff = fd.rate(pores=pn.pores('ymin'))

# %%
import napari
import porespy as ps
# im = ps.generators.overlapping_spheres(shape=[300, 300, 300], r=10, porosity=0.75)
napari.view_image(conc, rendering='attenuated_mip', attenuation=1.0, ndisplay=3, gamma=2, rotate=[10, 35, 0]);


# %% Compare tortuosity with dns
tau = ps.dns.tortuosity(im, axis=1, return_im=True)
plt.figure(2)
plt.imshow((tau.image/im).T)
print(tau.tortuosity)
