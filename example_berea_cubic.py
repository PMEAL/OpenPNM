import openpnm as op
import numpy as np
import scipy as sp
import openpnm.models as mods
import matplotlib.pyplot as plt


proj = op.materials.BereaCubic(shape=[15, 15, 15])
net = proj.network
geo = proj[1]

plt.hist(x=geo['throat.size']*1e6, bins=25,
         weights=geo['throat.volume']*(1e3)**3, edgecolor='k')
plt.hist(x=geo['pore.size_z']*1e6, bins=25,
         weights=geo['pore.volume']*(1e3)**3, edgecolor='k')

hg = op.phases.Mercury(network=net)
air = op.phases.Air(network=net)
water = op.phases.Water(network=net)
hg['throat.contact_angle'] = 140
hg['throat.surface_tension'] = 0.48
water['throat.viscosity'] = 0.0001

phys_air = op.physics.GenericPhysics(network=net, phase=air, geometry=geo)
phys_hg = op.physics.GenericPhysics(network=net, phase=hg, geometry=geo)
phys_water = op.physics.GenericPhysics(network=net, phase=water, geometry=geo)

# %%Simulate capillary pressure curve
mod = op.models.physics.capillary_pressure.washburn_slit
phys_hg.add_model(propname='throat.entry_pressure', model=mod)

mip = op.algorithms.Porosimetry(network=net)
mip.setup(phase=hg)
mip.set_inlets(net.pores(['top', 'bottom']))
mip.run(points=25)
mip.plot_intrusion_curve()

# %%Calculating permeability
mod = op.models.physics.hydraulic_conductance.hagen_poiseuille_slit
phys_water.add_model(propname='throat.hydraulic_conductance', model=mod)

alg = op.algorithms.StokesFlow(network=net, phase=water)
BC1_pores = net.pores('pore.front')
alg.set_value_BC(values=202650, pores=BC1_pores)
BC2_pores = net.pores('pore.back')
alg.set_value_BC(values=101325, pores=BC2_pores)
alg.run()
Q = alg.rate(pores=net.pores('front'))

Lx, Ly, Lz = net.spacing
Nx, Ny, Nz = net.shape
A = (Ly*Lz)*(Ny*Nz)
L = Lx*Nx
mu = sp.mean(water['throat.viscosity'])
Kxx = Q*mu*L/(A*101325)
print("The permeability coefficient is:", Kxx)

# %%Calculating porosity
Vp = geo['pore.volume'][net.Ps]
Vt = geo['throat.volume'][net.Ts]
Vps = sp.sum(Vp)
Vts = sp.sum(Vt)
Vt = Vps + Vts
Vb = Nx*Ny*Nz*Lx*Ly*Lz
e = Vt/Vb
print("The porosity is:", e)

# %%Calculating Formation Factor F
mod = op.models.physics.electrical_conductance.slit
phys_water.add_model(propname='throat.electrical_conductance', model=mod)

Om = op.algorithms.OhmicConduction(network=net, phase=water)
BC1_pores = net.pores('pore.front')
Om.set_value_BC(values=20, pores=BC1_pores)
BC2_pores = net.pores('pore.back')
Om.set_value_BC(values=0, pores=BC2_pores)
Om.run()
I = Om.rate(pores=net.pores('front'))

A = (Ly*Lz)*Ny*Nz
L = Lx*Nx
delta_V = 20
Rnet = delta_V*A/(I*L)
F = Rnet
print("The formation factor is:", F)
