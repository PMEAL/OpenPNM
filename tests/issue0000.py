import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)
import OpenPNM

#Place error causing lines below this line
#-----------------------------------------------------------------------------
pn = OpenPNM.Network.Cubic(division=[10,10,10],lattice_spacing=1,name='cubic1')