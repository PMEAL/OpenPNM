import OpenPNM as op
print('-----> Using OpenPNM version: '+op.__version__)

# =============================================================================
'''Build Topological Network'''
# =============================================================================
pn = op.Network.Cubic(shape=[50, 50, 1], spacing=0.0001, name='net')

# =============================================================================
'''Build Geometry'''
# =============================================================================
geom = op.Geometry.Toray090(network=pn, pores=pn.Ps, throats=pn.Ts)
geom.models['pore.seed']['seed'] = 0
geom['pore.seed'] = geom.models['pore.seed'].regenerate()
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

for i in range(0,20):
    qs = IP.queues
    for q in qs:
        IP.run2(queue=q, volume = 1e-12)
    IP.qregen()

plt.imshow(pn.asarray(IP['pore.invaded'])[:,:,0], interpolation='none')