# %% Package imports
import numpy as np
import openpnm as op
import matplotlib.pyplot as plt
import openpnm.models.geometry.diffusive_size_factors as gd
np.random.seed(10)

# Network
Nx = 100
shape = [Nx, Nx, 1]
spacing = 1e-2
net = op.network.Cubic(shape=shape, spacing=spacing)

# 2D square geometry (for comparison to comsol)
geo = op.geometry.GenericGeometry(network=net, pores=net.Ps, throats=net.Ts)
geo.add_model(propname='pore.diameter', model=op.models.misc.constant, value=spacing)
geo.add_model(propname='throat.length', model=op.models.misc.constant, value=1e-15)
geo.add_model(propname='throat.diameter', model=op.models.misc.constant, value=spacing)
geo["pore.cross_sectional_area"] = spacing
geo["throat.cross_sectional_area"] = spacing
geo["pore.volume"] = spacing**2
mod = gd.squares_and_rectangles
geo.add_model(propname='throat.diffusive_size_factors', model=mod)

# phase and physics
air = op.phase.Air(network=net)
phys = op.physics.GenericPhysics(network=net, phase=air, geometry=geo)

# Make diffusivity a linear function of temperature
air['pore.temperature'] = 300
air.add_model(propname='pore.diffusivity',
              model=op.models.misc.linear,
              m=1.860793056e-06,
              b=-0.0005375624384,
              prop='pore.temperature')
phys.add_model(propname='throat.diffusive_conductance',
               model=op.models.physics.diffusive_conductance.generic_diffusive)

# add constant thermal conductivity
air.remove_model(propname='pore.thermal_conductivity')
air["pore.thermal_conductivity"] = 0.0262

# add thermal conducance model to physics
phys.add_model(propname='throat.thermal_conductance',
               model=op.models.physics.thermal_conductance.generic_thermal)

tfd_settings = {
    "conductance": "throat.diffusive_conductance",
    "quantity": "pore.concentration",
}
tfc_settings = {
    "conductance": "throat.thermal_conductance",
    "quantity": "pore.temperature",
    "pore_volume": "pore.heat_capacity",
}

pardiso = op.solvers.PardisoSpsolve()
rk45 = op.integrators.ScipyRK45(verbose=True)

# Define Algorithms
# First algorithm, transient Fourier conduction
tfc = op.algorithms.TransientReactiveTransport(network=net, phase=air)
geo['pore.heat_capacity'] = geo['pore.volume'] * 1.0035 * 1000 * 1.225
tfc.settings._update(tfc_settings)
tfc.set_value_BC(net.pores("left"), 400)

# Second algorithm, transient Fickian diffusion
tfd = op.algorithms.TransientReactiveTransport(network=net, phase=air)
tfd.settings._update(tfd_settings)
tfd.set_value_BC(net.pores("left"), 100)

# Add variable props to algs
tfd.set_variable_props('pore.temperature')
tfc.set_variable_props('pore.concentration')

# handle initial conditions
T0 = np.ones(tfc.Np) * 300
T0[net.pores('left')] = 400
c0 = np.ones(tfd.Np)*50
c0[net.pores('left')] = 100
y0 = np.hstack((T0, c0))  # ICs must include boundary condition

# Integrator parameters
t_initial = 0
t_final = 50
t_step = 5
n_steps = int((t_final - t_initial)/t_step) + 1
t = np.linspace(t_initial, t_final, n_steps)
tspan = [t_initial, t_final]

# define transient multiphysics solver
tmp = op.contrib.TransientMultiPhysics(algorithms=[tfc, tfd], network=net)

# run
sol = tmp.run(y0, tspan, saveat=t)

# plot
T = sol[0:net.Np, -1]
C = sol[net.Np:, -1]

not_BC_pores = net.pores("left", mode='nor')
C_avg = sol[net.Np:, 0:][not_BC_pores].mean(axis=0)

# Load COMSOL data from .npz file
time = np.linspace(0, 2000, 401)
data = np.load('tmp_comsol_data.npz')
C_comsol = data['a']

# Plot concentration and temperature profiles
plt.figure(1)
fig, ax = plt.subplots(ncols=2)
im_1 = ax[0].imshow(T.reshape((Nx, Nx)))
im_2 = ax[1].imshow(C.reshape((Nx, Nx)))
fig.colorbar(im_1, ax=ax[0], fraction=0.046, pad=0.04)
fig.colorbar(im_2, ax=ax[1], fraction=0.046, pad=0.04)
ax[0].title.set_text('Temperature (K)')
ax[1].title.set_text('Concentration (mol m-3)')
ax[0].axis('off')
ax[1].axis('off')
im_1.set_clim(300, 400)
im_2.set_clim(0, 100)
fig.suptitle('OpenPNM', y=0.85)

# Plot average concentration of COMSOL compared to OpenPNM
plt.figure(2)
fig1, ax1 = plt.subplots()
ax1.plot(t, C_avg, time, C_comsol)
ax1.legend(('OpenPNM', 'COMSOL'))
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Concentration mol/m3')
ax1.set_title('Average Concentration')
