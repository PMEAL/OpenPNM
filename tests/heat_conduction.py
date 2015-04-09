import OpenPNM
import scipy as sp
import matplotlib.pylab as plt
from OpenPNM.Geometry import models as gm
from OpenPNM.Physics import models as pm

#==============================================================================
'''Build Topological Network'''
#==============================================================================
divs = [10,50]
Lc   = 0.1  # cm
pn = OpenPNM.Network.Cubic(name='net', shape= divs, spacing = Lc)
pn.add_boundaries()
#Trim z-direction boundaries
Ps = pn.pores(['top_boundary','bottom_boundary'])
pn.trim(pores=Ps)

#==============================================================================
'''Build Geometry'''
#==============================================================================
Ps = pn.pores('internal')
Ts = pn.throats()
geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=Ps,throats=Ts,name='geom')
geom['pore.area']     = Lc**2
geom['pore.diameter'] = Lc
geom['throat.length'] = 1e-25
geom['throat.area']   = Lc**2

Ps = pn.pores('boundary')
boun = OpenPNM.Geometry.GenericGeometry(network=pn,pores=Ps,name='boundary')
boun['pore.area']     = Lc**2
boun['pore.diameter'] =  1e-25

#==============================================================================
'''Build Material'''
#==============================================================================
Cu = OpenPNM.Phases.GenericPhase(network=pn,name='copper')
Cu['pore.thermal_conductivity'] = 1.0  # W/m.K

#==============================================================================
'''Build Physics'''
#==============================================================================
Ps = pn.pores()
Ts = pn.throats()
phys = OpenPNM.Physics.GenericPhysics(network=pn,phase=Cu,pores=Ps,throats=Ts)
phys.add_model(propname='throat.thermal_conductance',
               model=pm.thermal_conductance.series_resistors)

#==============================================================================
'''Run Algorithms'''
#==============================================================================
Fourier_alg = OpenPNM.Algorithms.FourierConduction(network=pn,phase=Cu)
inlets = pn.pores('back_boundary')
outlets = pn.pores(['front_boundary','left_boundary','right_boundary'])
T_out = 50  # Kelvin
T_in = 30*sp.sin(sp.pi*pn['pore.coords'][inlets,1]/5)+50
Fourier_alg.set_boundary_conditions(bctype='Dirichlet',bcvalue=T_in,pores=inlets)
Fourier_alg.set_boundary_conditions(bctype='Dirichlet',bcvalue=T_out,pores=outlets)
Fourier_alg.run()

Fourier_alg.return_results()

Cu['pore.analytical_temp'] = 30*sp.sinh(sp.pi*pn['pore.coords'][:,0]/5)/sp.sinh(sp.pi/5)*sp.sin(sp.pi*pn['pore.coords'][:,1]/5) + 50
a = Cu['pore.temperature'][pn.pores('geom')]
b = Cu['pore.analytical_temp'][pn.pores('geom')]
a = sp.reshape(a,(divs[0],divs[1]))
b = sp.reshape(b,(divs[0],divs[1]))
plt.subplot(3,1,1)
plt.imshow(b,interpolation='none')
plt.colorbar()
plt.subplot(3,1,2)
plt.imshow(a,interpolation='none')
plt.colorbar()
plt.subplot(3,1,3)
plt.imshow(a-b,interpolation='none')
plt.colorbar()
