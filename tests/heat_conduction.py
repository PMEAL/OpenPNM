import OpenPNM
import scipy as sp
import matplotlib.pylab as plt
from OpenPNM.Geometry import models as gm
from OpenPNM.Physics import models as pm

#==============================================================================
'''Build Topological Network'''
#==============================================================================
divs = [1,50,10]
Lc   = 0.1  # cm
pn = OpenPNM.Network.Cubic(name='net', shape= divs, spacing = Lc, loglevel=20)
pn.add_boundaries()
Ps = pn.pores(['front_boundary','back_boundary'])
pn.trim(pores=Ps)

#==============================================================================
'''Build Geometry'''
#==============================================================================
Ps = pn.pores('boundary',mode='difference')
Ts = pn.throats()
geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=Ps,throats=Ts,name='geom',loglevel=20)
geom['pore.area']     = Lc**2
geom['pore.diameter'] = Lc
geom['throat.length'] = 1e-25
geom['throat.area']   = Lc**2

Ps = pn.pores('boundary')
boun = OpenPNM.Geometry.GenericGeometry(network=pn,pores=Ps,name='boundary',loglevel=20)
boun['pore.area']     = Lc**2
boun['pore.diameter'] =  1e-25
boun['throat.area']   = Lc**2
boun['throat.length'] = 1e-25

#==============================================================================
'''Build Material'''
#==============================================================================
Cu = OpenPNM.Phases.GenericPhase(network=pn,name='copper',loglevel=20)
Cu['pore.thermal_conductivity'] = sp.nan
Cu['pore.thermal_conductivity'] = 1  # W/m.K

#==============================================================================
'''Build Physics'''
#==============================================================================
Ps = pn.pores()
Ts = pn.throats()
phys = OpenPNM.Physics.GenericPhysics(network=pn,phase=Cu,pores=Ps,throats=Ts,loglevel=10)
phys.add_model(propname='throat.thermal_conductance',
               model=pm.thermal_conductance.series_resistors)

phys.regenerate()  # Update the conductance values
# To prevent insulating the network by very low conductance
phys['throat.thermal_conductance'][pn.find_neighbor_throats(pn.pores('boundary'))] = sp.average(phys['throat.thermal_conductance'][pn.find_neighbor_throats(pn.pores('boundary',mode='difference'))])
#==============================================================================
'''Run Algorithms'''
#==============================================================================
Fourier_alg = OpenPNM.Algorithms.FourierConduction(network=pn,phase=Cu,loglevel=10)
inlets = pn.pores('top_boundary')
outlets = pn.pores(['bottom_boundary','left_boundary','right_boundary'])
T_out = 50  # Kelvin
T_in = 30*sp.sin(sp.pi*pn['pore.coords'][inlets,1]/5)+50
Fourier_alg.set_boundary_conditions(bctype='Dirichlet',bcvalue=T_in,pores=inlets)
Fourier_alg.set_boundary_conditions(bctype='Dirichlet',bcvalue=T_out,pores=outlets)
Fourier_alg.run()

Fourier_alg.update_results()

import OpenPNM.Utilities.IO as io
io.VTK.save(network=pn,phases=Cu)

Cu['pore.analytical_temp'] = 30*sp.sinh(sp.pi*pn['pore.coords'][:,2]/5)/sp.sinh(sp.pi/5)*sp.sin(sp.pi*pn['pore.coords'][:,1]/5) + 50
a = Cu['pore.temperature'][pn.pores('geom')]
b = Cu['pore.analytical_temp'][pn.pores('geom')]
a = sp.reshape(a,(divs[2],divs[1]))
b = sp.reshape(b,(divs[2],divs[1]))
plt.imshow(a-b,interpolation='none')
plt.colorbar()

