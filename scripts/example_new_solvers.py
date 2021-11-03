# %% Initializations
import numpy as np
import openpnm as op
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from openpnm.utils import tic, toc

np.random.seed(10)
ws = op.Workspace()
ws.clear()
ws.settings["loglevel"] = 40

# %% Set the stage up for algorithms
Nx = 10
shape = [Nx, Nx, Nx]
spacing = 1/Nx
net = op.network.Cubic(shape=shape, spacing=spacing)
air = op.phases.Air(network=net)
air["throat.diffusive_conductance"] = np.random.rand(net.Nt)

# Add reactions
mod = op.models.physics.generic_source_term.standard_kinetics
k, n = -1e-1, 2
air["pore.k"], air["pore.n"] = k, n
air.add_model(
    X="pore.concentration",
    propname="pore.rxn",
    model=mod, prefactor="pore.k", exponent="pore.n",
    regen_mode="deferred"
)

alg_settings = {
    "conductance": "throat.diffusive_conductance",
    "quantity": "pore.concentration"
}
pardiso = op.solvers.PardisoSpsolve()
rk45 = op.integrators.ScipyRK45(verbose=True)

# %% Try w/ GenericTransport
print(f"\n{'+'*65}", "\n--> Using GenericTransport", f"\n{'+'*65}\n", flush=True)
gt = op.algorithms.GenericTransport(network=net, phase=air)
gt.settings.update(alg_settings)
gt.set_value_BC(net.pores("left"), 1.0)
gt.set_value_BC(net.pores("right"), 0.4)

gt.run(solver=pardiso)

print(gt["pore.concentration"].mean())

# %% Try w/ ReactiveTransport
print(f"\n{'+'*65}", "\n--> Using ReactiveTransport", f"\n{'+'*65}\n", flush=True)
rt = op.algorithms.ReactiveTransport(network=net, phase=air)
rt.settings.update(alg_settings)

rt.set_value_BC(net.pores("left"), 1.0)
rt.set_value_BC(net.pores("right"), 0.4)
rt.set_source(propname="pore.rxn", pores=net.pores("surface", mode="not"))

rt.run(solver=pardiso)

print(rt["pore.concentration"].mean())

# %% Try w/ FickianDiffusion
print(f"\n{'+'*65}", "\n--> Using FickianDiffusion", f"\n{'+'*65}\n", flush=True)

fd = op.algorithms.FickianDiffusion(network=net, phase=air)

fd.set_value_BC(net.pores("left"), 1.0)
fd.set_value_BC(net.pores("right"), 0.4)

fd.run(solver=pardiso)

print(fd["pore.concentration"].mean())

# %% Playing around w/ Transient stuff
print(f"\n{'+'*65}", "\n--> Using TransientReactiveTransport", f"\n{'+'*65}\n", flush=True)

# air["throat.diffusive_conductance"] = np.ones(net.Nt, dtype=float) * spacing
air["throat.diffusive_conductance"] = np.random.rand(net.Nt) * spacing
V = net["pore.volume"] = spacing**3

trt = op.algorithms.TransientReactiveTransport(network=net, phase=air)
trt.settings.update(alg_settings)
trt.set_value_BC(net.pores("left"), 1.0)
trt.set_value_BC(net.pores("right"), 0.4)
c0 = np.zeros(trt.Np)

tspan = [0, 0.4]
tout = np.linspace(tspan[0], tspan[1], 11)
dt = (tspan[1] - tspan[0]) / 200

# Adding source term
trt.set_source("pore.rxn", pores=net.pores("surface", mode="not"))
# air.remove_model("pore.rxn")

tic()
sol = trt.run(x0=c0, tspan=tspan, integrator=rk45, saveat=tout)
toc()

# %% Hardcoded solution using solve_ivp
# trt._build_A()
# trt._build_b()
# trt._apply_BCs()
# trt._apply_sources()
# A = trt.A.tocsc()
# b = trt.b

# from scipy.integrate import solve_ivp

# def ode_func(t, y, A, b, V):
#     return (-A.dot(y) + b) / V

# f = lambda t, y: ode_func(t, y, A, b, V)

# tic()
# # TODO: implicit methods: pass lband/uband for drastically better performance
# # TODO: lband/uband = 1 is probably wrong!
# sol = solve_ivp(fun=f, t_span=tspan, y0=c0, method="RK23", rtol=1e-5)
# toc()

# %% TransientReactiveTransport (legacy)
from openpnm.algorithms.legacy import TransientReactiveTransport
trt_legacy = TransientReactiveTransport(network=net, phase=air)
dt_avg = np.diff(sol.t).mean()
t_settings = {
    "t_initial": tspan[0],
    "t_final": tspan[1],
    "t_scheme": "implicit",
    "t_step": dt,
    "t_output": tout.tolist()
}
trt_legacy.settings.update({**alg_settings, **t_settings})
trt_legacy.set_value_BC(net.pores("left"), 1.0)
trt_legacy.set_value_BC(net.pores("right"), 0.4)
trt_legacy.set_IC(0.0)

# Add source term
trt_legacy.set_source("pore.rxn", pores=net.pores("surface", mode="not"))

tic()
trt_legacy.run(solver=pardiso)
toc()

# %% Compare new implementation vs legacy
c_new = sol.reshape((*shape, sol.t.size))

c_legacy = np.zeros((trt_legacy.Np, tout.size))
for i, t in enumerate(tout):
    t_fmt = trt_legacy._nbr_to_str(t)
    c_legacy[..., i] = trt_legacy[f"pore.concentration@{t_fmt}"]
c_legacy = c_legacy.reshape((*shape, tout.size))

fig, ax = plt.subplots(ncols=2)
ax[0].imshow(c_new[..., -1].mean(axis=2))
ax[0].set_title("new implementation")
ax[1].set_title("legacy")
h = ax[1].imshow(c_legacy[..., -1].mean(axis=2))

fig, ax = plt.subplots()
ax.plot(sol.t, c_new.mean(axis=(0, 1, 2)), "g-", label="new implementation")
ax.plot(tout, c_legacy.mean(axis=(0, 1, 2)), "rx", label="legacy")
ax.legend()

# Ensure correctness
np.testing.assert_allclose(c_new[..., -1], c_legacy[..., -1], rtol=1e-3)

# %% Continue integration
# # Integrate from 0 to 0.2 at once
# sol1 = trt.run(x0=c0, tspan=[0, 0.2], solver=rk45)

# # Integrate from 0 to 0.1 and then from 0.1 to 0.2
# sol2 = trt.run(x0=c0, tspan=[0, 0.1], solver=rk45)
# sol3 = trt.run(x0=sol2(0.1), tspan=[0.1, 0.2], solver=rk45)

# # Ensure correctness
# np.testing.assert_allclose(sol1(0.2), sol3(0.2), rtol=1e-5)
