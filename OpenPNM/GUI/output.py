import os
import OpenPNM

def VTK(net, filename='new_network.vtp'):
  filename = str(filename)
  OpenPNM.Visualization.VTK().write( **locals() )
  