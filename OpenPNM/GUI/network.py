from __future__ import absolute_import
import numpy as np
from .. import GEN

def Cubic(divisions=[10,10,10], lattice_spacing=.0005, loglevel=20, btype=[0,0,0], C=0.9, D=2.09e-9):
  gn = GEN.Cubic(divisions=divisions, lattice_spacing=lattice_spacing, loglevel=loglevel, btype=btype)
  pn = gn.generate()
  setattr(pn,"Total_Conc",C)
  setattr(pn,"Diff_Coefficient",D)
  setattr(pn,"divisions",divisions)
  setattr(pn,"lattice_spacing",lattice_spacing)
  pn.throat_properties['Pc_entry'] = -4*0.072*np.cos(np.radians(105))/pn.throat_properties['diameter']
  return {'network': pn}