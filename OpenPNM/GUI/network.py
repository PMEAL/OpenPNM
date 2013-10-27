from __future__ import absolute_import
import os
import numpy as np
import OpenPNM

def Cubic(domain_size=[1,1,1], divisions=[10,10,10], stats_pores=None, stats_throats=None, btype=[0,0,0]):
  pn = OpenPNM.Geometry.Cubic().generate(**locals())
  return {'net': pn}