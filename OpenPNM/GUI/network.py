from __future__ import absolute_import
from .. import GEN

def Cubic(divisions=[10,10,10],lattice_spacing=.0005,loglevel=20,btype = [0,0,0]):
  gn = GEN.Cubic(divisions=divisions, lattice_spacing=lattice_spacing, loglevel=loglevel, btype=btype)
  return {'network':gn.generate()}