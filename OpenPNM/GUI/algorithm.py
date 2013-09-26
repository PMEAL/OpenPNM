from __future__ import absolute_import
from .. import Algorithms
import numpy as np

def OrdinaryPercolation(network, npts=100, inv_faces=[1], P=0.04):
  OP = Algorithm.OrdinaryPercolation(network, npts=npts, inv_faces=inv_faces)
  OP.run()
  network.set_pore_property('INVADED PORES',np.multiply(np.ones(network.get_num_pores()),network.pore_properties['Pc_invaded']<P))
  Algorithms.FickianDiffusion(network, Alg='OP', Pressure=[P])
  return {'network': network}