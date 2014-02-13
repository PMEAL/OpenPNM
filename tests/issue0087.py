import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM
import sp

#Place error causing lines below this line
#-----------------------------------------------------------------------------
pn = OpenPNM.Network.TestNet()
loc = [0,1]
pn.set_pore_info(label='test',locations=loc) #this fails
pn.set_pore_info(label='test',locations=sp.array(loc)) #this works