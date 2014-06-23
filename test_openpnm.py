#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Tests for OpenPNM """

from __future__ import print_function, division, absolute_import
import pytest

import OpenPNM
import numpy as np
import scipy as sp

def test_template_generator():
    R = np.array([[[0,0],[0,1]]])
    pn = OpenPNM.Network.Template()
    pn.generate(R)
    pn.trim(R>0.5)
    assert len(pn.get_pore_data(prop='coords'))==3
    assert len(pn.get_throat_data(prop='conns'))==2

#def test_linear_solver():
#    # fix cube dimensions?
#    pn = OpenPNM.Network.Cubic()
#    pn.generate(add_boundaries=False)
#
#    x,y,z = pn.get_pore_data(prop='coords').T
#    left = x==x.min()
#    right = x==x.max()
#    ics = 2*left + 1*right
#
#    sol = OpenPNM.Shortcuts.solve_linear(pn, ics)
#    assert( np.allclose(sol, np.array([2,1.5,1]*9)) )

#def test_rectilinear_integrity_after_prune():
#    R = np.random.rand(50,40,30)
#    # some simple visual pruning, for comparison
#    M = np.where(R > R.mean(), R, 0)
#    # the convoluted graph way
#    pn = OpenPNM.Network.Template()
#    pn.generate(R)
#    pn.trim(R<=R.mean())
#    O = pn.asarray()
#    # what it would look like normally
#    assert np.allclose(M, O)
    
def test_graph_queries():
    pn = OpenPNM.Network.TestNet()
    pn.clone(pores=[0])
    pn.num_pores()
    assert pn.find_neighbor_pores(pores=pn.pores()[-1]) == [0]
    assert sp.all(pn.find_neighbor_pores(pores=[0]) == [  1,   5,  25, 125])
    pn.stitch(heads=[124],tails=[125])
    assert pn.num_throats() == 302
    pn.extend(pore_coords=[[1,1,1]])
    assert pn.find_neighbor_pores(pores=126) == []
    assert pn.find_neighbor_throats(pores=126) == []    
    
def test_transport():
    pn = OpenPNM.Network.TestNet()
    geom = OpenPNM.Geometry.Stick_and_Ball(network=pn)
    geom.set_locations(pores=pn.pores('all'), throats=pn.throats('all'))
    geom.add_property(prop='pore_seed',model='random',seed=2)
    geom.regenerate()
    air = OpenPNM.Fluids.Air(network=pn)
    air.regenerate()
    phys = OpenPNM.Physics.GenericPhysics(network=pn,fluid=air,geometry=geom)
    phys.add_method(prop='diffusive_conductance',model='bulk_diffusion')
    phys.regenerate()
    fick = OpenPNM.Algorithms.FickianDiffusion(network=pn)
    fick.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.6,pores=pn.pores('top'))
    fick.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.4,pores=pn.pores('bottom'))
    fick.run(active_fluid=air)
    fick.update()
    data = sp.round_(air['pore.mole_fraction'][[0,25,50,75,100]],6)
    assert sp.sum(sp.in1d(data, [ 0.4, 0.461799, 0.519957, 0.570009, 0.6])) == 5
    

if __name__ == '__main__':
    pytest.main([__file__])