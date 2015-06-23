import OpenPNM
import pytest
from OpenPNM.Utilities import topology
ctrl = OpenPNM.Base.Controller()
topo = topology()


def test_subdivide():
    pn = OpenPNM.Network.Cubic(shape=[5, 5, 5],
                               spacing=0.001,
                               name='micro_net')
    pn['pore.micro'] = True
    nano_pores = [2, 13, 14, 15]
    pn.subdivide(pores=nano_pores, shape=[4, 4, 4], labels='nano')
    assert pn.Np == (125+4*64-4)
    assert pn.Nt == (300+(4*144)-16+15*16+16)
    ctrl.export(network=pn, filename='nano')


def test_clone_and_trim():
    ctrl.clear()
    pn = OpenPNM.Network.Cubic(shape=[5, 5, 5], name='net')
    geom = OpenPNM.Geometry.GenericGeometry(network=pn, name='geo1')
    geom.set_locations(pores=pn.Ps, throats=pn.Ts)
    assert sorted(list(ctrl.keys())) == ['geo1', 'net']
    pn2 = ctrl.clone_simulation(pn, name='clone')
    assert sorted(list(ctrl.keys())) == ['geo1', 'geo1_clone', 'net',
                                         'net_clone']
    topo.trim(network=pn2, pores=pn2.pores('top'))
