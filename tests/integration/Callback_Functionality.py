import numpy as np
import openpnm as op
import matplotlib.pyplot as plt


def test_callback_functionality():

    Nx = 101
    shape = [Nx, 1, 1]
    spacing = 1/Nx * 5
    net = op.network.Cubic(shape=shape, spacing=spacing)
    air = op.phase.Air(network=net)
    air["throat.diffusive_conductance"] = spacing * np.ones(net.Nt)
    net["pore.volume"] = spacing**3

    # Set up transient Fickian diffusion algorithm
    tfd = op.algorithms.TransientFickianDiffusion(network=net, phase=air)
    tfd.set_value_BC(net.pores("left"), 0)
    tfd.set_value_BC(net.pores("right"), 0)

    # Define a pulse signal
    def pulse(t, y):
        if 0 <= t <= 0.05:
            y[net.Np//2] = 1.0

    # Add the pulse signal to the algorithm as a callback
    tfd._set_callback(pulse)

    # Solve the transient algorithm
    c0 = np.zeros(tfd.Np)
    tspan = [0, 0.4]
    tfd.run(x0=c0, tspan=tspan)

    # Plot c vs. time
    tout = np.linspace(tspan[0], tspan[1], 10)
    fig, ax = plt.subplots()
    for i, t in enumerate(tout):
        ax.plot(tfd.soln['pore.concentration'](t), label=f"{t:.2f} (s)")
    ax.legend()
    ax.set_title("Dissipation of a pulse signal , c(x=0,L) = 0")
    ax.set_xlabel("distance (m)")
    ax.set_ylabel("concentration (mol/m$^3$)")
