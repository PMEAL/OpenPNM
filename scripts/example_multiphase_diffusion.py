r"""
Example: Multiphase diffusion with heterogeneous reaction

2D network, consists of air and water. Air occupies the middle of the
network and is surrounded by two film-like regions of water at the top
and the bottom. The top and the bottom faces of the network is assumed
to be coated w/ a catalyst, and are therefore reactive.

The diffusing species diffuses through air, which is followed by mass
partitioning at the two air-water interfaces, then continues diffusing
through the water, and finally reacting at the two reacting plates at the
top and the bottom of the network.

"""
import openpnm as op
from openpnm.network import Cubic
from openpnm.phases import Air, Water, MultiPhase
from openpnm.geometry import StickAndBall
from openpnm.physics import Standard
from openpnm.models.misc import constant
from openpnm.models.physics.source_terms import linear
import matplotlib.pyplot as plt
import numpy as np
np.random.seed(10)

# Define network, geometry and constituent phases
shape = np.array([100, 100, 1])
net = Cubic(shape=shape)
geom = StickAndBall(network=net, pores=net.Ps, throats=net.Ts)
air = Air(network=net, name="air")
water = Water(network=net, name="water")
water["pore.diffusivity"] = air["pore.diffusivity"] * 0.05

# Define the regions to be occupied by the two phases (air and water)
x, y, z = net["pore.coords"].T
ps_water = net.Ps[(y >= 75) + (y <= 25)]
ps_air = np.setdiff1d(net.Ps, ps_water)
ts_water = net.find_neighbor_throats(pores=ps_water, mode="xnor")
ts_air = net.find_neighbor_throats(pores=ps_air, mode="xnor")
ts_interface = net.find_neighbor_throats(pores=ps_water, mode="xor")

# Define multiphase and set phase occupancy
mphase = MultiPhase(network=net, phases=[air, water], name="mphase")
mphase._set_automatic_throat_occupancy()
mphase.set_occupancy(air, pores=ps_air, throats=ts_air)
mphase.set_occupancy(water, pores=ps_water, throats=ts_water)

# Define physics
phys = Standard(network=net, phase=mphase, geometry=geom)
# Assign a partition coefficient (concentration ratio)
K_water_air = 0.5   # c @ water / c @ air
mphase.set_binary_partition_coef(phases=[water, air], model=constant, value=K_water_air)
# Replace the "default" ordinary_diffusion w/ multiphase_diffusion conductance model
mdiff = op.models.physics.diffusive_conductance.multiphase_diffusion
phys.add_model(propname="throat.diffusive_conductance", model=mdiff)

# Fickian diffusion
fd = op.algorithms.FickianDiffusion(network=net, phase=mphase)
# Set source term
phys["pore.A1"] = -1e-8 * geom["pore.area"]
phys["pore.A2"] = 0.0
phys.add_model(
    propname="pore.rxn",
    model=linear,
    X="pore.concentration",
    A1="pore.A1", A2="pore.A2",
    regen_mode="deferred"
)
rxn_pores = net.pores(["front", "back"])
net.set_label("rxn", pores=rxn_pores)
fd.set_source(propname="pore.rxn", pores=rxn_pores)

# Set BCs and run simulation
net.set_label("air", pores=ps_air)
left_air = net.pores(["left", "air"], mode="and")
right_air = net.pores(["right", "air"], mode="and")
fd.set_value_BC(pores=left_air, values=1.0)
fd.set_value_BC(pores=right_air, values=0.1)
fd.run()

# Post-processing
mphase.update(fd.results())
c = mphase["pore.concentration"]
c2d = np.rot90(c.reshape(shape).squeeze())
plt.imshow(c2d)
plt.colorbar()

op.io.XDMF.export_data(network=net, phases=mphase, filename="network")
