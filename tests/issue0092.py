import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM

#Place error causing lines below this line
#-----------------------------------------------------------------------------
pn = OpenPNM.Network.TestNet()
geom = OpenPNM.Geometry.Stick_and_Ball(network=pn, name='stick_and_ball', locations=pn.get_pore_indices())

# Add vector method
geom.add_method(prop='throat_vector', model='pore_to_pore')

# Generate geometry data
geom.regenerate()

print(pn.get_throat_data(prop='tvecs'))  # This is right, so getter is ok
print(pn.get_throat_data(prop='vector'))  # This is wrong, so setter is wrong
