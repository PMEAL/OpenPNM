from __future__ import absolute_import
import OpenPNM
import numpy as np

def OrdinaryPercolation(net, npts=100, inv_sites=[1]):
  OpenPNM.Algorithms.OrdinaryPercolation(**locals()).run()
  return {'net': net}