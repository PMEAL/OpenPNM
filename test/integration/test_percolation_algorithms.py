import scipy as sp
import OpenPNM
import pytest


def test_IP_old_approach():
    pn = OpenPNM.Network.Cubic(shape=[30, 30, 1], spacing=0.0001)
    geom = OpenPNM.Geometry.Toray090(network=pn)
    water = OpenPNM.Phases.Water(network=pn)
    phys = OpenPNM.Physics.Standard(network=pn, phase=water, geometry=geom)
    inlets = pn.pores('left')
    IP_1 = OpenPNM.Algorithms.InvasionPercolation(network=pn)
    IP_1.run(phase=water, inlets=inlets)
    a = ['pore.all', 'pore.invaded', 'pore.invasion_sequence', 'throat.all',
         'throat.entry_pressure', 'throat.invaded', 'throat.invasion_sequence',
         'throat.order', 'throat.sorted']
    assert sorted(list(IP_1.keys())) == a
    IP_1.return_results()
    assert 'throat.invasion_sequence' in water.keys()
    assert 'pore.invasion_sequence' in water.keys()


def test_IP_new_approach():
    pn = OpenPNM.Network.Cubic(shape=[30, 30, 1], spacing=0.0001)
    geom = OpenPNM.Geometry.Toray090(network=pn)
    water = OpenPNM.Phases.Water(network=pn)
    phys = OpenPNM.Physics.Standard(network=pn, phase=water, geometry=geom)
    inlets = pn.pores('left')
    IP_1 = OpenPNM.Algorithms.InvasionPercolation(network=pn)
    IP_1.setup(phase=water)
    IP_1.set_inlets(pores=inlets)
    IP_1.run()
    a = ['pore.all', 'pore.invaded', 'pore.invasion_sequence', 'throat.all',
         'throat.entry_pressure', 'throat.invaded', 'throat.invasion_sequence',
         'throat.order', 'throat.sorted']
    assert sorted(list(IP_1.keys())) == a
    IP_1.return_results()
    assert 'throat.invasion_sequence' in water.keys()
    assert 'pore.invasion_sequence' in water.keys()
