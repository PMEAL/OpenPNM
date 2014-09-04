import OpenPNM
import scipy as sp
from OpenPNM.Physics import models as pm

#==============================================================================
'''Build Topological Network'''
#==============================================================================
divs = [1,50,10]
Lc = 0.00004
pn = OpenPNM.Network.Cubic(name='net',shape = divs, spacing = Lc,loglevel=20)
pn.add_boundaries()

Ps = pn.pores(['front_boundary','back_boundary'])
pn.trim(pores=Ps)

#==============================================================================
'''Build Geometry'''
#==============================================================================
Ps = pn.pores('boundary',mode='difference')
Ts = pn.throats()
geom = OpenPNM.Geometry.Toray090(network=pn,pores=Ps,throats=Ts,name='GDL')

Ps = pn.pores('boundary')
Ts = pn.throats('boundary')
boun = OpenPNM.Geometry.Boundary(network=pn,pores=Ps,throats=Ts,name='boundary')

#==============================================================================
'''Build Material'''
#==============================================================================
air = OpenPNM.Phases.Air(network=pn,name='air')

#==============================================================================
'''Build Physics'''
#==============================================================================
Ps = pn.pores()
Ts = pn.throats()
phys = OpenPNM.Physics.GenericPhysics(network=pn,phase=air,pores=Ps,throats=Ts)
phys.add_model(propname='throat.hydraulic_conductance',
               model=pm.hydraulic_conductance.hagen_poiseuille)
phys.regenerate()  # Update the conductance values

#==============================================================================
'''Run Algorithms'''
#==============================================================================
Darcy = OpenPNM.Algorithms.StokesFlow(network=pn,phase=air)
inlets = pn.pores('bottom_boundary')
Ps = pn.pores('top_boundary')
outlets = Ps[pn['pore.coords'][Ps,1]<(divs[1]*Lc/2)]
P_out = 0  # Pa
Q_in = 0.6667*(Lc**2)*divs[1]*divs[0]  # m^3/s
Darcy.set_boundary_conditions(bctype='Neumann_group',bcvalue=-Q_in,pores=inlets)
Darcy.set_boundary_conditions(bctype='Dirichlet',bcvalue=P_out,pores=outlets)
Darcy.run()
Darcy.update_results()

Darcy2 = OpenPNM.Algorithms.StokesFlow(network=pn,phase=air)
inlets = pn.pores('bottom_boundary')
outlets = pn.pores('top_boundary')
P_out = 10  # Pa
P_in = 1000  # Pa
Darcy2.set_boundary_conditions(bctype='Dirichlet',bcvalue=P_in,pores=inlets)
Darcy2.set_boundary_conditions(bctype='Dirichlet',bcvalue=P_out,pores=outlets)

Darcy2.run()
Q = -Darcy2.rate(inlets)
K = Q*air['pore.viscosity'][0]*divs[2]*Lc/(divs[0]*divs[1]*Lc**2*(P_in-P_out))

Vp = sp.sum(pn['pore.volume']) + sp.sum(pn['throat.volume'])
Vb = sp.prod(divs)*Lc**3
e = Vp/Vb

import OpenPNM.Utilities.IO as io
io.VTK.save(network=pn,phases=[air])

