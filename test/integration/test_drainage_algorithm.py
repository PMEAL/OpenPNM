import scipy as sp
import matplotlib
import OpenPNM
import pytest

#def test_drainage_basic():
pn = OpenPNM.Network.Cubic(shape=[25, 25, 25], spacing=0.0001)
pn['pore.internal'] = True
pn['throat.internal'] = True
pn.add_boundary_pores(pores=pn.pores('top'),
                      offset=[0, 0, 0.0001],
                      apply_label='boundary_top')
pn.add_boundary_pores(pores=pn.pores('bottom'),
                      offset=[0, 0, -0.0001],
                      apply_label='boundary_bottom')
geom = OpenPNM.Geometry.Toray090(network=pn,
                                 pores=pn.pores('internal'),
                                 throats=pn.throats('internal'))
boun = OpenPNM.Geometry.Boundary(network=pn,
                                 pores=pn.pores('internal', mode='not'),
                                 throats=pn.throats('internal', mode='not'))
water = OpenPNM.Phases.Water(network=pn)
phys = OpenPNM.Physics.Standard(network=pn, phase=water, geometry=geom)

drainage = OpenPNM.Algorithms.Drainage(network=pn)
drainage.setup(invading_phase=water)
drainage.set_inlets(pores=pn.pores('boundary_top'))
drainage.run()
data = drainage.get_drainage_data()
assert sp.amin(data['nonwetting_phase_saturation']) == 0.0
assert sp.amax(data['nonwetting_phase_saturation']) == 1.0
