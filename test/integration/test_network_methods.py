import scipy as sp
import OpenPNM
import pytest

def test_find_connected_pores():
    pn = OpenPNM.Network.Cubic(shape=(10,10,10))
    a = pn.find_connected_pores(throats=[0,1])
    assert sp.all(a.flatten() == [0, 1, 1, 2])
    a = pn.find_connected_pores(throats=[0,1], flatten=True)
    assert sp.all(a == [0, 1, 2])
    Tind = sp.zeros((pn.Nt,),dtype=bool)
    Tind[[0,1]] = True
    a = pn.find_connected_pores(throats=Tind, flatten=True)
    assert sp.all(a == [0, 1, 2])
    a = pn.find_connected_pores(throats=[], flatten=True)
    assert sp.shape(a) == (0,2)

def test_find_neighbor_pores():
    pn = OpenPNM.Network.Cubic(shape=(10,10,10))
    a = pn.find_neighbor_pores(pores=[])
    assert sp.size(a) == 0
    Pind = sp.zeros((pn.Np,), dtype=bool)
    Pind[[0,1]] = True
    a = pn.find_neighbor_pores(pores=Pind)
    assert sp.all(a == [2, 10, 11, 100, 101])
    a = pn.find_neighbor_pores(pores=[0, 2], mode='union')
    assert sp.all(a == [1, 3, 10,  12, 100, 102])
    a = pn.find_neighbor_pores(pores=[0, 2], mode='intersection')
    assert sp.all(a == [1])
    a = pn.find_neighbor_pores(pores=[0, 2], mode='not_intersection')
    assert sp.all(a == [3, 10, 12, 100, 102])
    a = pn.find_neighbor_pores(pores=[0, 2], mode='union', excl_self=False)
    assert sp.all(a == [ 0, 1, 2, 3, 10, 12, 100, 102])
    a = pn.find_neighbor_pores(pores=[0, 2], mode='intersection', excl_self=False)
    assert sp.all(a == [1])
    a = pn.find_neighbor_pores(pores=[0, 2], mode='not_intersection', excl_self=False)
    assert sp.all(a == [0, 2, 3, 10, 12, 100, 102])
