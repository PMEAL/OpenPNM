import openpnm as op
import scipy as sp
import scipy.sparse.csgraph as spgr

ws = op.core.Workspace()
ws.settings['local_data'] = True

########################## NETWORK ##########################
sp.random.seed(0)
pn = op.network.Cubic(shape=[120, 30, 1], spacing=2e-5, name='pn11')
#pn.add_boundary_pores()

########################## GEOMETRIES #########################
geom = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)

########################## PHASES ##########################
water = op.phases.Water(network=pn)

########################## PHYSICS ##########################
phys_water = op.physics.GenericPhysics(network=pn, phase=water, geometry=geom)

water['throat.viscosity'] = water['pore.viscosity'][0]
mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys_water.add_model(propname='throat.hydraulic_conductance',
                     model=mod, viscosity='throat.viscosity')

water['pore.diffusivity'] = 2.14e-5
water['pore.molar_density'] = 1
geom['pore.area'] = sp.pi*(geom['pore.diameter']**2)/4.0
mod2 = op.models.physics.diffusive_conductance.bulk_diffusion
phys_water.add_model(propname='throat.diffusive_conductance',
                     model=mod2, diffusivity='pore.diffusivity')

phys_water.regenerate_models()

inlet = pn.pores('back') # pore inlet
outlet = pn.pores('front') # pore outlet

inlet2 = pn.pores('left') # pore inlet2
outlet2 = pn.pores('right') # pore outlet2

########################## ALGORITHMS ##########################


alg1=op.algorithms.StokesFlow(network=pn,phase=water)
alg1.setup()
alg1.set_BC(pores=inlet2, bctype='dirichlet', bcvalues=1.5*101325)
alg1.set_BC(pores=outlet2, bctype='dirichlet', bcvalues=101325)
alg1['pore.pressure'] = 101325
alg1.run()
water['pore.pressure'] = alg1['pore.pressure']


#mod3 = op.models.physics.advective_diffusive_conductance.advection_diffusion
#phys_water.add_model(propname='throat.advective_diffusive_conductance',
#                     model=mod3, diffusivity='pore.diffusivity',
#                     hydraulic_conductance='throat.hydraulic_conductance',
#                     pressure='pore.pressure')
#phys_water.regenerate_models()
#am = pn.create_adjacency_matrix(weights=-water['throat.advective_diffusive_conductance'], \
#                                fmt='coo')
#A3 = spgr.laplacian(am)
#A3 = A3.todense()


alg2=op.algorithms.FickianDiffusion(network=pn,phase=water)
alg2.setup()
alg2.set_BC(pores=inlet, bctype='dirichlet', bcvalues=1)
alg2.set_BC(pores=outlet, bctype='dirichlet', bcvalues=0)
alg2.run()



alg3=op.algorithms.AdvectionDiffusion(network=pn,phase=water)
alg3.setup()
alg3.set_BC(pores=inlet, bctype='dirichlet', bcvalues=1)
alg3.set_BC(pores=outlet, bctype='dirichlet', bcvalues=0)
alg3.run()



#CONNS = pn['throat.conns'][[water['throat.all']]]
#TS = sp.array(pn.find_connecting_throat(CONNS[:,0],CONNS[:,1]))
#
#Q = water['throat.hydraulic_conductance'][TS]*\
#                   sp.diff(water['pore.pressure'][CONNS],axis=1)[:,0]
#
#QP1 = sp.where(Q>0, Q, 0)
#QN1 = sp.where(Q<0, Q, 0)
#
#QP2 = -QN1
#QN2 = -QP1
#
##QP2 = sp.where(-Q>0, -Q, 0)
##QN2 = sp.where(-Q<0, -Q, 0)
#
#g1 = sp.empty((water.Nt, 2), dtype=float)
#g1[:,0] = QN1-water['throat.diffusive_conductance']
#g1[:,1] = QN2-water['throat.diffusive_conductance']
#am1 = pn.create_adjacency_matrix(weights=-g1, fmt='coo')
#A1 = spgr.laplacian(am1)
#
#g2 = sp.empty((water.Nt, 2), dtype=float)
#g2[:,0] = -QP1-water['throat.diffusive_conductance']
#g2[:,1] = -QP2-water['throat.diffusive_conductance']
#am2 = pn.create_adjacency_matrix(weights=-g2, fmt='coo')
#A2 = spgr.laplacian(am2)
#
#A2.setdiag(A1.diagonal())
#
#A2 = A2.todense()

#################################################


#A3 = sp.zeros((pn.Np,pn.Np))
#b = sp.zeros(shape=(pn.Np, ), dtype=float)
#
#nt = pn.find_neighbor_throats(pores=pn['pore.all'], \
#     flatten=False, mode='not_intersection') # pore neighbor throats
#
#HC = water['throat.hydraulic_conductance'] # throat hydraulic conductance
#p = water['pore.pressure'] # pore pressure
#D = water['throat.diffusive_conductance'] # throat diffusive conductance
#
#for i in range (pn.Np):
#    q = HC[nt[i]]*(p[pn.find_neighbor_pores(i)]-p[i])
#    qP = sp.where(q>0, q, 0)
#    qN = sp.where(q<0, q, 0)
#    
#    A3[i,i] = sp.sum( qN - D[nt[i]] ) # mehmani
#
#    j1 = pn['throat.conns'][nt[i]]
#    j2 = sp.reshape(j1,sp.size(j1))
#    j = j2[j2!=i]
#
#    A3[i,j] = -( -qP - D[nt[i]] ) # mehmani
#
#
#
#
#
#
#
#
#
#A3[inlet,:] = 0
#A3[inlet,inlet] = 1
#b[inlet] = 1
#
#A3[outlet, :] = 0
#A3[outlet, outlet] = 1
#b[outlet] = 0
#
#A3 = sp.sparse.csr_matrix(A3)
#x3 = sp.sparse.linalg.spsolve(A3, b)
#water['pore.conc_mehmani3'] = x3

#A2[inlet,:] = 0
#A2[inlet,inlet] = 1
#
#A2[outlet, :] = 0
#A2[outlet, outlet] = 1
#
#A2 = sp.sparse.csr_matrix(A2)
#x2 = sp.sparse.linalg.spsolve(A2, b)
#water['pore.conc_mehmani2'] = x2


water['pore.concD'] = alg2['pore.mole_fraction']
water['pore.concAD'] = alg3['pore.mole_fraction']

ws.export_data(network=pn, phases=water, filename='vis_test')








