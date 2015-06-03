import OpenPNM
import scipy as sp
import matplotlib.pyplot as plt
print('-----> Using OpenPNM version: '+OpenPNM.__version__)

ctrl = OpenPNM.Base.Controller()
#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic(shape=[25,1,100],spacing=0.001,name='net')
pn.add_boundaries()

#==============================================================================
'''Build Geometry'''
#==============================================================================
Ps = pn.pores('boundary',mode='not')
Ts = pn.find_neighbor_throats(pores=Ps,mode='intersection',flatten=True)
geom = OpenPNM.Geometry.Toray090(network=pn,pores=Ps,throats=Ts)
geom.models['pore.seed']['seed'] = 0
geom.models['pore.seed']['regen_mode'] = 'normal'
geom.regenerate()

Ps = pn.pores('boundary')
Ts = pn.find_neighbor_throats(pores=Ps,mode='not_intersection')
boun = OpenPNM.Geometry.Boundary(network=pn,pores=Ps,throats=Ts)

#==============================================================================
'''Build Phases'''
#==============================================================================
water = OpenPNM.Phases.Water(network=pn,name='water')
air = OpenPNM.Phases.Air(network=pn,name='air')
air['pore.surface_tension'] = 0.072
air['pore.contact_angle'] = 0

#==============================================================================
'''Build Physics'''
#==============================================================================
Ps = pn.pores()
Ts = pn.throats()
phys_water = OpenPNM.Physics.Standard(network=pn,phase=water,pores=Ps,throats=Ts)
#Add some additional models to phys_water
phys_water.models.add(model=OpenPNM.Physics.models.capillary_pressure.static_pressure,
                      propname='pore.static_pressure',
                      regen_mode='deferred')
phys_air = OpenPNM.Physics.Standard(network=pn,phase=air,pores=Ps,throats=Ts)
phys_air.models['throat.capillary_pressure'] = phys_water.models['throat.capillary_pressure'].copy()

#==============================================================================
'''Begin Simulations'''
#==============================================================================
'''Perform Invasion Percolation using Version 2'''
#------------------------------------------------------------------------------
IP = OpenPNM.Algorithms.InvasionPercolation2(network=pn)
inlets = pn.pores('top_boundary')
filled = False
count = 0
#while not filled:
while count < 20:
    count += 1
    print('Loop number:', count)
    IP.setup(phase=air, p_inlets=inlets)
    filled = IP.run(nsteps=50)
    water['throat.occupancy'] = False
    water['throat.occupancy'][IP['throat.invaded'] == -1] = True
    phys_water.models.regenerate('pore.static_pressure')
    P12 = pn['throat.conns']
    temp = sp.amax(phys_water['pore.static_pressure'][P12],axis=1)
    phys_water['throat.capillary_pressure'] += temp
    inlets = (IP['pore.invaded'] >= 0)

Ps = pn.pores('internal')
plt.matshow(pn.asarray(phys_water['pore.static_pressure'][Ps])[:,0,:].T,
                       interpolation='none',
                       origin='lower')

plt.matshow(pn.asarray(IP['pore.invaded'][Ps])[:,0,:].T,
                       interpolation='none',
                       origin='lower')























