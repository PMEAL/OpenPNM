import scipy as sp
import OpenPNM
import pytest


def test_find_connected_pores():
    pn = OpenPNM.Network.Cubic(shape=(10, 10, 10))
    a = pn.find_connected_pores(throats=[0, 1])
    assert sp.all(a.flatten() == [0, 1, 1, 2])
    a = pn.find_connected_pores(throats=[0, 1], flatten=True)
    assert sp.all(a == [0, 1, 2])
    Tind = sp.zeros((pn.Nt,), dtype=bool)
    Tind[[0, 1]] = True
    a = pn.find_connected_pores(throats=Tind, flatten=True)
    assert sp.all(a == [0, 1, 2])
    a = pn.find_connected_pores(throats=[], flatten=True)
    assert sp.shape(a) == (0, )


def test_find_neighbor_pores():
    pn = OpenPNM.Network.Cubic(shape=(10, 10, 10))
    a = pn.find_neighbor_pores(pores=[])
    assert sp.size(a) == 0
    Pind = sp.zeros((pn.Np,), dtype=bool)
    Pind[[0, 1]] = True
    a = pn.find_neighbor_pores(pores=Pind)
    assert sp.all(a == [2, 10, 11, 100, 101])
    a = pn.find_neighbor_pores(pores=[0, 2], mode='union')
    assert sp.all(a == [1, 3, 10,  12, 100, 102])
    a = pn.find_neighbor_pores(pores=[0, 2], mode='intersection')
    assert sp.all(a == [1])
    a = pn.find_neighbor_pores(pores=[0, 2], mode='not_intersection')
    assert sp.all(a == [3, 10, 12, 100, 102])
    a = pn.find_neighbor_pores(pores=[0, 2],
                               mode='union',
                               excl_self=False)
    assert sp.all(a == [0, 1, 2, 3, 10, 12, 100, 102])
    a = pn.find_neighbor_pores(pores=[0, 2],
                               mode='intersection',
                               excl_self=False)
    assert sp.all(a == [1])
    a = pn.find_neighbor_pores(pores=[0, 2],
                               mode='not_intersection',
                               excl_self=False)
    assert sp.all(a == [0, 2, 3, 10, 12, 100, 102])


def test_find_neighbor_throats():
    pn = OpenPNM.Network.Cubic(shape=(10, 10, 10))
    a = pn.find_neighbor_throats(pores=[])
    assert sp.size(a) == 0
    Pind = sp.zeros((pn.Np,), dtype=bool)
    Pind[[0, 1]] = True
    a = pn.find_neighbor_throats(pores=Pind)
    assert sp.all(a == [0, 1, 900, 901, 1800, 1801])
    a = pn.find_neighbor_throats(pores=[0, 2], mode='union')
    assert sp.all(a == [0, 1, 2, 900, 902, 1800, 1802])
    a = pn.find_neighbor_throats(pores=[0, 2], mode='intersection')
    assert sp.size(a) == 0
    a = pn.find_neighbor_throats(pores=[0, 2], mode='not_intersection')
    assert sp.all(a == [0, 1, 2, 900, 902, 1800, 1802])


def test_num_neighbors():
    pn = OpenPNM.Network.Cubic(shape=(10, 10, 10))
    a = pn.num_neighbors(pores=[])
    assert sp.size(a) == 0
    Pind = sp.zeros((pn.Np,), dtype=bool)
    Pind[0] = True
    a = pn.num_neighbors(pores=Pind)
    assert a == 3
    a = pn.num_neighbors(pores=[0, 2], flatten=True)
    assert a == 6
    assert isinstance(a, int)
    a = pn.num_neighbors(pores=[0, 2], flatten=False)
    assert sp.all(a == [3, 4])
    a = pn.num_neighbors(pores=0, flatten=False)
    assert sp.all(a == [3])
    assert isinstance(a, sp.ndarray)


def test_find_interface_throats():
    pn = OpenPNM.Network.Cubic(shape=(3, 3, 3))
    pn['pore.domain1'] = False
    pn['pore.domain2'] = False
    pn['pore.domain3'] = False
    pn['pore.domain1'][[0, 1, 2]] = True
    pn['pore.domain2'][[5, 6, 7]] = True
    pn['pore.domain3'][18:26] = True
    a = pn.find_interface_throats(labels=['domain1', 'domain2'])
    assert a == [20]
    a = pn.find_interface_throats(labels=['domain1', 'domain3'])
    assert sp.size(a) == 0
