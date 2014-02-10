import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM

#Place error causing lines below this line
#-----------------------------------------------------------------------------
pn = OpenPNM.Network.TestNet()
loc = pn.get_pore_indices()
geo = OpenPNM.Geometry.GenericGeometry(name='geo_test',network=pn,locations=loc)
geo.add_method(prop='pore_seed',model='constant',value=0.123)
seeds = pn.get_pore_data(subdomain='geo_test',prop='pore_seed')
type(seeds) #this is empty
pn._pore_data.keys() #pore_seed has not been added