import OpenPNM
import scipy as sp
from OpenPNM.Geometry import models as gm
from OpenPNM.Physics import models as pm

#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic(name='net',loglevel=20)
divs = [1,50,10]
Lc   = 0.1  # cm
pn.generate(divisions=divs,lattice_spacing=[Lc],add_boundaries=True)

Ps = pn.pores(['front_face','back_face'])
pn.trim(pores=Ps)

#==============================================================================
'''Build Geometry'''
#==============================================================================
Ps = pn.pores('boundary',mode='difference')
Ts = pn.throats()
geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=Ps,throats=Ts,name='geom')
geom['pore.area']     = Lc**2
geom['pore.diameter'] = Lc
geom['throat.length'] = 1e-25
geom['throat.area']   = Lc**2

Ps = pn.pores('boundary')
boun = OpenPNM.Geometry.GenericGeometry(network=pn,pores=Ps,name='boundary')
boun['pore.area']     = Lc**2
boun['pore.diameter'] = 1e-25
boun['throat.area']   = Lc**2
boun['throat.length'] = 1e-25

#==============================================================================
'''Build Material'''
#==============================================================================
Cu = OpenPNM.Fluids.GenericFluid(network=pn,name='copper')
Cu['pore.thermal_conductivity'] = sp.nan
Cu['pore.thermal_conductivity'] = 1  # W/m.K

#==============================================================================
'''Build Physics'''
#==============================================================================
Ps = pn.pores()
Ts = pn.throats()
phys = OpenPNM.Physics.GenericPhysics(network=pn,fluid=Cu,pores=Ps,throats=Ts)
phys.add_model(propname='throat.thermal_conductance',
               model=pm.thermal_conductance.series_resistors)


#==============================================================================
'''Run Algorithms'''
#==============================================================================
Keff = OpenPNM.Algorithms.FourierConduction(network=pn,fluid=Cu)
inlets = pn.pores('top_face')
outlets = pn.pores(['bottom_face','left_face','right_face'])
T_out = 50  # Kelvin
T_in = 30*sp.sin(sp.pi*pn['pore.coords'][inlets,1]/5)+50
Keff.set_boundary_conditions(bctype='Dirichlet',bcvalue=T_in,pores=inlets)
Keff.set_boundary_conditions(bctype='Dirichlet',bcvalue=T_out,pores=outlets)

phys.regenerate()  # Update the conductance values
Keff.setup(fluid=Cu)
Keff.run()

Keff.update_results()

#vis = OpenPNM.Visualization.VTK()
#vis.write(network=pn,fluids=[Cu])

Cu['pore.analytical_temp'] = 30*sp.sinh(sp.pi*pn['pore.coords'][:,2]/5)/sp.sinh(sp.pi/5)*sp.sin(sp.pi*pn['pore.coords'][:,1]/5) + 50
a = Cu['pore.temperature'][pn.pores('geom')]
b = Cu['pore.analytical_temp'][pn.pores('geom')]
a = sp.reshape(a,(divs[2],divs[1]))
b = sp.reshape(b,(divs[2],divs[1]))
plt.imshow(a-b,interpolation='none')
plt.colorbar()

