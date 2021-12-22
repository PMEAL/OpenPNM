import scipy as sp
import numpy as np
import openpnm as op
import pandas as pd
import matplotlib.pyplot as plt
import openpnm.models.geometry.diffusive_size_factors as gd
import timeit

np.random.seed(10)

# %% Test scipy's solve_ivp
# FOR EXAMPLE: solve dy / dt = y^2 + y

def fun(t, y):
    
    return y

sol = sp.integrate.solve_ivp(fun, t_span=(0, 10), y0=np.array([5]))

# plot
t = sol.t
y = np.reshape(sol.y, newshape=len(sol.t))
plt.plot(t, y)

# %% Set up for the solvers
Nx = 100
shape = [Nx, Nx, 1]
spacing = 1e-2 #1/Nx # adjust 
net = op.network.Cubic(shape=shape, spacing=spacing)
# geo = op.geometry.CirclesAndRectangles(network=net, pores=net.Ps, throats=net.Ts)
# 2d square geometry
geo = op.geometry.GenericGeometry(network=net, pores=net.Ps, throats=net.Ts)
geo.add_model(propname='pore.diameter',
               model=op.models.misc.constant,
               value=spacing)
geo.add_model(propname='throat.length',
               model=op.models.misc.constant,
               value=1e-15)
geo.add_model(propname='throat.diameter',
               model=op.models.misc.constant,
               value=spacing)
Ndim = op.topotools.dimensionality(network=net).sum()
geo.add_model(propname='pore.cross_sectional_area',
               model=op.models.misc.constant,
               value=spacing**(Ndim-1))
geo.add_model(propname='throat.cross_sectional_area',
               model=op.models.misc.constant,
               value=spacing**(Ndim-1))
geo.add_model(propname='pore.volume',
               model=op.models.misc.constant,
               value=spacing**Ndim)
mod = gd.squares_and_rectangles if Ndim == 2 else gd.cubes_and_cuboids
geo.add_model(propname='throat.diffusive_size_factors',
               model=mod,
               pore_diameter="pore.diameter",
               throat_diameter="throat.diameter")

# Diff = op.metrics.EffectiveDiffusivity(network=net)
# # Diff.settings._update({
# # 'area': (10*1e-4)**(Ndim-1)})
# #'length':9*1e-4,})
# D_eff = Diff.run()
# print(f"Effective diffusivity in x direction is: {D_eff}")

air = op.phases.Air(network=net)
phys = op.physics.GenericPhysics(network=net, phase=air, geometry=geo)

# make diffusivity a function of temperature - ALREADY IS!!
air['pore.temperature'] = 300
air.add_model(propname='pore.diffusivity',
              model=op.models.misc.linear, 
              m=1.860793056e-06,
              b=-0.0005375624384,
              prop='pore.temperature')
'''
air.add_model(propname='pore.diffusivity',
              model=op.models.misc.constant,
              value=2.067547840000001e-04)
'''
phys.add_model(propname='throat.diffusive_conductance', 
               model=op.models.physics.diffusive_conductance.generic_diffusive)

air.remove_model(propname='pore.thermal_conductivity')
air.add_model(propname='pore.thermal_conductivity',
              model=op.models.misc.constant,
              value=0.0262,
              regen_mode='constant')

phys.add_model(propname='throat.thermal_conductance',
               model=op.models.physics.thermal_conductance.generic_thermal) 

tfd_settings = {
    "conductance": "throat.diffusive_conductance",
    "quantity": "pore.concentration",
    "cache_A": False,
    "cache_b": False  
}

tfc_settings = {
    "conductance": "throat.thermal_conductance",
    "quantity": "pore.temperature",
    "pore_volume": "pore.heat_capacity",
    "cache_A": False,
    "cache_b": False  
}

pardiso = op.solvers.PardisoSpsolve()
rk45 = op.integrators.ScipyRK45(verbose=True)

# %% Define Algorithms
# First algorithm, transient fourier conduction
tfc = op.algorithms.TransientReactiveTransport(network=net, phase=air)
geo['pore.heat_capacity'] = geo['pore.volume'] * 1.0035 * 1000 * 1.225
tfc.settings._update(tfc_settings)
tfc.set_value_BC(net.pores("left"), 400)
T0 = np.ones(tfc.Np) * 300
T0[net.pores('left')] = 400

# Second algorithm, transient fickian diffusion
tfd = op.algorithms.TransientReactiveTransport(network=net, phase=air)
# tfd.settings['variable_props'] = 'pore.temperature'
tfd.settings._update(tfd_settings)
tfd.set_value_BC(net.pores("left"), 100)
# tfd.set_value_BC(net.pores("right"), 100)
c0 = np.ones(tfd.Np)*50
c0[net.pores('left')] = 100

# time
t_initial = 0
t_final = 50
t_step = 5
n_steps = int((t_final - t_initial)/t_step) + 1
t = np.linspace(t_initial, t_final, n_steps)

# %% For Loop Solution
# solve multiphysics system assuming temperature change is small over t_step
not_BC_pores = net.pores("left", mode='nor')  
t_prev = 0
C_1_avg = [c0[not_BC_pores].mean()]
for ti in t:
    if ti == t_initial:
        continue
    print('time:', ti, "s")
    tspan = [t_prev, ti]
    t_prev = ti
    tout = ti
    # Solve for temperature first... add while loop
    # temperature dependent thermal conductivity
    sol_1 = tfc.run(x0=T0, tspan=tspan, integrator=rk45, saveat=tout)
    air.regenerate_models() # update diffusivuty
    phys.regenerate_models() # update diffusive conductance because tfd does not have iterative props
    sol_2 = tfd.run(x0=c0, tspan=tspan, integrator=rk45, saveat=tout)
    # update initial coniditions
    T0 = sol_1[:, 1]
    c0 = sol_2[:, 1]
    C_1_avg.append(c0[not_BC_pores].mean())

# note: dual coupling, time step needs to be small, hamed benchamrk solution
T_1 = sol_1[:, -1]
C_1 = sol_2[:, -1]

# %% Build RHS manually
# add variable props to algs
tfd.set_variable_props('pore.temperature')
tfc.set_variable_props('pore.concentration')

def _build_rhs(algs):
    
    def ode_func(t, y):
        
        # initialize rhs
        rhs = []
        for i, alg in enumerate(algs):
            # get x from y, assume alg.Np is same for all algs
            x = y[i*tfc.Np:(i+1)*tfc.Np]
            # store x onto algorithm
            alg.x = x
            # build A and b
            alg._update_iterative_props()
            alg._build_A()
            alg._build_b()
            alg._apply_BCs()
            alg._apply_sources()
            A = alg.A.tocsc()
            b = alg.b
            V = geo[alg.settings['pore_volume']]
            rhs_alg = np.hstack(-A.dot(x) + b)/V
            rhs = np.hstack((rhs, rhs_alg))

        return rhs
       
    return ode_func

# call solve_ivp and pass function that builds rhs as dydt
algs = [tfc, tfd]
rhs = _build_rhs(algs)
T0 = np.ones(tfc.Np) * 300
T0[net.pores('left')] = 400
c0 = np.ones(tfd.Np)*50
c0[net.pores('left')] = 100

y0 = np.hstack((T0, c0)) # ICs must include boundary condition
tspan = [0, t_final]
rtol = 1e-5
# sol = sp.integrate.solve_ivp(fun=rhs, t_span=tspan, y0=y0, method="RK45", t_eval=t, rtol=rtol)
tmp = op.algorithms.TransientMultiPhysics(algorithms=[tfc, tfd], network=net)
sol = tmp.run(y0, tspan, saveat=t)

# plot
# T_2 = sol.y[0:net.Np, -1]
# C_2 = sol.y[net.Np:, -1]
# C_2_avg = sol.y[net.Np:, 0:][not_BC_pores].mean(axis=0)
T_2 = sol.T.transpose()[0:net.Np, -1]
C_2 = sol.T.transpose()[net.Np:, -1]
C_2_avg = sol.T.transpose()[net.Np:, 0:][not_BC_pores].mean(axis=0)

# %% Pure Fickian Diffusion and Pure Fourier Conduction
air.remove_model(propname='pore.diffusivity')
air.add_model(propname='pore.diffusivity',
              model=op.models.misc.constant,
              value=2.067547840000001e-04)
phys.regenerate_models()
tspan = [t_initial, t_final]
sol_tfd = tfd.run(x0=c0, tspan=tspan, integrator=rk45, saveat=t)
sol_tfc = tfc.run(x0=T0, tspan=tspan, integrator=rk45, saveat=t)

C_pure = sol_tfd.T[:, 0:][:, not_BC_pores].mean(axis=1)
T_pure = sol_tfc.T[:, 0:][:, not_BC_pores].mean(axis=1)

# %% Import COMSOL data and make plots
io = r'C:\Users\mmcka\OneDrive - University of Waterloo\UW files\code\MultiphysicsSolver\heat&mass.xlsx'
df = pd.read_excel(io, sheet_name=['Coupled','Diffusion', 'Conduction'])
time = np.asarray(df['Coupled'].iloc[0:, 0])
C_3_avg = np.asarray(df['Coupled'].iloc[0:, 1])
C_pure_comsol = np.asarray(df['Diffusion'].iloc[0:, 1])
T_pure_comsol = np.asarray(df['Conduction'].iloc[0:, 1])

# plot average concentrations comparing for loop, build rhs, and comsol
fig1, ax1 = plt.subplots()
# ax1.plot(t, C_1_avg, time, C_3_avg)
# ax1.legend(('For Loop', 'COMSOL'))
ax1.plot(t, C_1_avg, t, C_2_avg, time, C_3_avg)
ax1.legend(('For Loop', 'Build RHS', 'COMSOL'))
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Concentration mol/m3')
ax1.set_title('Average Concentration')

# plot average concentration from pure fickian diffusion
fig2, ax2 = plt.subplots()
ax2.plot(t, C_pure, time, C_pure_comsol)
ax2.legend(('OpenPNM', 'COMSOL'))
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Average Concentration (mol/m3)')
ax2.set_title('Only Diffusion')

# plot average concentration from pure fourier conduction
fig3, ax3 = plt.subplots()
ax3.plot(t, T_pure, time, T_pure_comsol)
ax3.legend(('OpenPNM', 'COMSOL'))
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Average Temperature (mol/m3)')
ax3.set_title('Only Conduction')

plt.figure(1)
fig, ax = plt.subplots(ncols=2)
im_1 = ax[0].imshow(T_1.reshape((Nx,Nx)))
im_2 = ax[1].imshow(C_1.reshape((Nx,Nx)))
fig.colorbar(im_1, ax=ax[0], fraction=0.046, pad=0.04)
fig.colorbar(im_2, ax=ax[1], fraction=0.046, pad=0.04)
ax[0].title.set_text('Temperature (K)')
ax[1].title.set_text('Concentration (mol m-3)')
plt.axis('off')
im_1.set_clim(300, 400)
im_2.set_clim(0, 100)
fig.suptitle('For Loop', y=0.85)

plt.figure(2)
fig, ax = plt.subplots(ncols=2)
im_1 = ax[0].imshow(T_2.reshape((Nx,Nx)))
im_2 = ax[1].imshow(C_2.reshape((Nx,Nx)))
fig.colorbar(im_1, ax=ax[0], fraction=0.046, pad=0.04)
fig.colorbar(im_2, ax=ax[1], fraction=0.046, pad=0.04)
ax[0].title.set_text('Temperature (K)')
ax[1].title.set_text('Concentration (mol m-3)')
plt.axis('off')
im_1.set_clim(300, 400)
im_2.set_clim(0, 100)
fig.suptitle('Build RHS', y=0.85)

# Error calculation
T_error = np.abs(T_2 - T_1) / T_2 * 100
C_error = np.abs(C_2 - C_1) / C_2 * 100
print(T_error.max(), C_error.max())

plt.figure(3)
fig, ax = plt.subplots(ncols=2)
im_1 = ax[0].imshow(T_error.reshape((Nx,Nx)))
im_2 = ax[1].imshow(C_error.reshape((Nx,Nx)))
fig.colorbar(im_1, ax=ax[0], fraction=0.046, pad=0.04)
fig.colorbar(im_2, ax=ax[1], fraction=0.046, pad=0.04)
ax[0].title.set_text('Temperature (K)')
ax[1].title.set_text('Concentration (mol m-3)')
plt.axis('off')
# im_1.set_clim(0, 1)
# im_2.set_clim(0, 1)
fig.suptitle('Error', y=0.85)

# %% Export Geometry
# op.io.COMSOL.export_data(network=net, filename='multiphysics_solver_2d')
# comparison for 2D square