import scipy as sp
import OpenPNM
import pytest
ctrl = OpenPNM.Base.Controller()


def test_IP_old_approach():
    ctrl.clear()
    pn = OpenPNM.Network.Cubic(shape=[30, 30, 1], spacing=0.01)
    geom = OpenPNM.Geometry.Toray090(network=pn, pores=pn.Ps, throats=pn.Ts)
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
    ctrl.clear()


def test_IP_new_approach():
    ctrl.clear()
    pn = OpenPNM.Network.Cubic(shape=[30, 30, 1], spacing=0.0001)
    geom = OpenPNM.Geometry.Toray090(network=pn, pores=pn.Ps, throats=pn.Ts)
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
    ctrl.clear()


def test_OP_old_approach():
    ctrl.clear()
    pn = OpenPNM.Network.Cubic(shape=[30, 30, 1], spacing=0.0001)
    geom = OpenPNM.Geometry.Toray090(network=pn, pores=pn.Ps, throats=pn.Ts)
    water = OpenPNM.Phases.Water(network=pn)
    phys = OpenPNM.Physics.Standard(network=pn, phase=water, geometry=geom)
    OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=pn,
                                                  invading_phase=water)
    Ps = pn.pores(labels=['left'])
    OP_1.run(inlets=Ps)
    OP_1.return_results(Pc=7000)
    a = ['pore.all', 'pore.inlets', 'pore.inv_Pc', 'pore.inv_sat',
         'pore.inv_seq', 'throat.all', 'throat.entry_pressure',
         'throat.inv_Pc', 'throat.inv_sat', 'throat.inv_seq']
    assert sorted(list(OP_1.keys())) == a
    ctrl.clear()


def test_OP_new_approach():
    ctrl.clear()
    pn = OpenPNM.Network.Cubic(shape=[30, 30, 1], spacing=0.0001)
    geom = OpenPNM.Geometry.Toray090(network=pn, pores=pn.Ps, throats=pn.Ts)
    water = OpenPNM.Phases.Water(network=pn)
    phys = OpenPNM.Physics.Standard(network=pn, phase=water, geometry=geom)
    inlets = pn.pores('left')
    OP = OpenPNM.Algorithms.OrdinaryPercolation(network=pn)
    OP.setup(invading_phase=water)
    OP.set_inlets(pores=inlets)
    OP.run(npts=25)
    a = ['pore.all', 'pore.inlets', 'pore.inv_Pc', 'pore.inv_sat',
         'pore.inv_seq', 'throat.all', 'throat.entry_pressure',
         'throat.inv_Pc', 'throat.inv_sat', 'throat.inv_seq']
    assert sorted(list(OP.keys())) == a
    ctrl.clear()
