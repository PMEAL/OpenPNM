from __future__ import absolute_import
import OpenPNM
import numpy as np

def OrdinaryPercolation(network, npts=100, inv_faces=[1], P=0.04):
  OP = OpenPNM.Algorithms.OrdinaryPercolation(network, npts=npts, inv_faces=inv_faces)
  OP.run()
  network.set_pore_property('INVADED PORES',np.multiply(np.ones(network.get_num_pores()),network.pore_properties['Pc_invaded']<P))
  OpenPNM.Algorithms.FickianDiffusion(network, Alg='OP', Pressure=[P])
  return {'network': network}
  
def InvasionPercolation(network, end_condition='breakthrough',timing='ON',report=20):
    mask = network.pore_properties['inlets'] == 1   
    inlets = network.pore_properties['numbering'][mask]
    mask = network.pore_properties['outlets'] == 1   
    outlets = network.pore_properties['numbering'][mask]
    IP = OpenPNM.Algorithms.InvasionPercolation(network, inlets=inlets, outlets=outlets, end_condition=end_condition, timing=timing, report=report)
    IP.run()
    return {'network': network}