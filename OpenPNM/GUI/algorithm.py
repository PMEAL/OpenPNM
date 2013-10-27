from __future__ import absolute_import
import OpenPNM
import scipy as sp

def OrdinaryPercolation(net, invading_fluid, defending_fluid,
                        npts=50, AL=True):
  invading_fluid.set_pair(defending_fluid)
  invading_fluid.regenerate()
  defending_fluid.regenerate()

  OP_1 = OpenPNM.Algorithms.OrdinaryPercolation()
  OpenPNM.Physics.CapillaryPressure.Washburn(net,invading_fluid)
  a = net.pore_properties['type']==1
  OP_1.run(network=net,
           invading_fluid=invading_fluid,
           defending_fluid=defending_fluid,
           inv_sites=a,npts=npts,AL=AL)

  return {'net':net}

def InvasionPercolation(net, invading_fluid, defending_fluid,
                        in_pore_type=3,
                        out_pore_type=4,
                        quarter_lim=0.1,
                        npts=50, AL=True):
  invading_fluid.set_pair(defending_fluid)
  invading_fluid.regenerate()
  defending_fluid.regenerate()

  IP_1 = OpenPNM.Algorithms.InvasionPercolation()
  OpenPNM.Physics.CapillaryPressure.Washburn(net, defending_fluid)
  face = net.pore_properties['type']==in_pore_type
  quarter = sp.rand(net.get_num_pores(),)<quarter_lim
  inlets = net.pore_properties['numbering'][face&quarter]
  outlets = net.pore_properties['numbering'][net.pore_properties['type']==out_pore_type]
  IP_1.run(net,
           invading_fluid=invading_fluid,
           defending_fluid=defending_fluid,
           inlets=inlets,outlets=outlets)
  return {'net':net}

def DiffusionSimulation(net, invading_fluid, defending_fluid, Pc=4000):
  OpenPNM.Physics.MultiPhase.update_occupancy_OP(invading_fluid, Pc)
  OpenPNM.Physics.MultiPhase.effective_occupancy(net, defending_fluid)
  OpenPNM.Physics.MassTransport.DiffusiveConductance(net, defending_fluid)

  Fickian_alg = OpenPNM.Algorithms.FickianDiffusion()

  BCtypes = sp.zeros(net.get_num_pores())
  BCvalues = sp.zeros(net.get_num_pores())

  BCtypes[net.pore_properties['type']==1] = 1
  BCtypes[net.pore_properties['type']==6] = 1
  BCvalues[net.pore_properties['type']==1] = 8e-2
  BCvalues[net.pore_properties['type']==6] = 8e-1

  Fickian_alg.set_boundary_conditions(types=BCtypes,values=BCvalues)

  Fickian_alg.run(net, active_fluid=defending_fluid)
