import OpenPNM as op
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation

print('-----> Using OpenPNM version: '+op.__version__)
ctrl = op.Base.Controller()
ctrl.loglevel = 30

# =============================================================================
'''Build Topological Network'''
# =============================================================================
pn = op.Network.Cubic(shape=[30, 30, 1], spacing=1, name='net')

# =============================================================================
'''Build Geometry'''
# =============================================================================
geom = op.Geometry.Toray090(network=pn, pores=pn.Ps, throats=pn.Ts)
geom.models['pore.seed']['seed'] = 0
geom['pore.seed'] = geom.models['pore.seed'].regenerate()
del geom.models['pore.diameter']
geom['pore.diameter'] = geom['pore.seed'].copy()
geom.models.add(propname='throat.diameter',
                model=op.Geometry.models.throat_misc.neighbor,
                mode='min')
geom.models.regenerate()

# =============================================================================
'''Build Phases'''
# =============================================================================
air = op.Phases.Air(network=pn, name='air')
water = op.Phases.Water(network=pn, name='water')

# =============================================================================
'''Build Physics'''
# =============================================================================
Ps = pn.pores()
Ts = pn.throats()
phys_water = op.Physics.Standard(network=pn, phase=water, pores=Ps, throats=Ts)
phys_air = op.Physics.Standard(network=pn, phase=air, pores=Ps, throats=Ts)

# =============================================================================
'''Begin Simulations'''
# =============================================================================

# -----------------------------------------------------------------------------
'''Perform Invasion Percolation'''
# -----------------------------------------------------------------------------
IP = op.Algorithms.InvasionPercolationDrying(network=pn)
IP.setup(invading_phase=air, defending_phase=water)
IP.set_inlets(pores=pn.pores('front'))

D = op.Algorithms.FickianDiffusion(network=pn, phase=air)
Ps_in = pn.pores('front')

ims = []
fig = plt.figure()
i = 0
while sp.sum(IP['pore.invaded'] == -1) > 0:
    # Calculate the diffusion of vapor given the existing invasion pattern
    Ps_out = pn.toindices(IP['pore.invaded'] == -1)
    Ps_in = pn.pores('front')
    D.set_boundary_conditions(pores=Ps_out,
                              bcvalue=0.5,
                              bctype='Dirichlet',
                              mode='overwrite')
    D.set_boundary_conditions(pores=Ps_in,
                              bcvalue=0.1,
                              bctype='Dirichlet',
                              mode='merge')
    D.run()

    # Update the invasion pattern based on vapor loss from each cluster
    qs = IP.queues
    for qnum in range(len(qs)):
        print(str(qnum)+' : ', end='')
        r = sp.absolute(D.rate(pores=IP['pore.queue_number'] == qnum))
        IP.run2(queue=qs[qnum], volume=r*6000)
    IP.qregen()

    # Create image and store for making an animation
    im = sp.copy(D['pore.air_mole_fraction'])
    im[Ps_out] = sp.nan
    im = pn.asarray(im)[:, :, 0]
    ims.append((plt.imshow(im, interpolation='none'),))
    # Purge diffusion object, so it can be restarted (fix this)
    i += 1

im_ani = animation.ArtistAnimation(fig,
                                   ims,
                                   interval=50,
                                   repeat_delay=3000,
                                   blit=True)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# im_ani.save('test2.mp4', writer=writer)
