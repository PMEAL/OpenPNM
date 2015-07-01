import OpenPNM as op
import matplotlib.pyplot as plt
plt.ion()

print('-----> Using OpenPNM version: '+op.__version__)
ctrl = op.Base.Controller()
ctrl.loglevel = 30

# =============================================================================
'''Build Topological Network'''
# =============================================================================
pn = op.Network.Cubic(shape=[50, 50, 1], spacing=1, name='net')

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

for i in range(0, 402):
    qs = IP.queues
    for i in range(len(qs)):
        q = qs[i]
        IP.run2(queue=q, volume=0.5)
    IP.qregen()

plt.imshow(pn.asarray(IP['pore.queue_number'])[:, :, 0], interpolation='none')
