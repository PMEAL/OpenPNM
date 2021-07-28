import numpy as np
import openpnm as op
np.random.seed(10)

net = op.network.Cubic(shape=[10, 10, 10])
air = op.phases.Air(network=net)
air["throat.diffusive_conductance"] = np.random.rand(net.Nt)

alg_settings = {
    "conductance": "throat.diffusive_conductance",
    "quantity": "pore.concentration"
}
solver = op.solvers.ScipySpsolve()

# %% Try w/ GenericTransport
print("--Using GenericTransport")
gt = op.algorithms.GenericTransport(network=net, phase=air)
gt.settings.update(alg_settings)
gt.set_value_BC(net.pores("left"), 1.0)
gt.set_value_BC(net.pores("right"), 0.4)

gt.run(solver=solver)

print(gt["pore.concentration"].mean())

# %% Try w/ ReactiveTransport
print("--Using ReactiveTransport")
rt = op.algorithms.ReactiveTransport(network=net, phase=air)
rt.settings.update(alg_settings)

rt.set_value_BC(net.pores("left"), 1.0)
rt.set_value_BC(net.pores("right"), 0.4)

rt.run(solver=solver)

print(rt["pore.concentration"].mean())

# %% Try w/ FickianDiffusion
print("--Using FickianDiffusion")

fd = op.algorithms.FickianDiffusion(network=net, phase=air)

fd.set_value_BC(net.pores("left"), 1.0)
fd.set_value_BC(net.pores("right"), 0.4)

fd.run(solver=solver)

print(fd["pore.concentration"].mean())
