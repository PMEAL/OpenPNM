#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Tests for OpenPNM """

from __future__ import print_function, division, absolute_import
import pytest

import OpenPNM
import numpy as np

def test_template_generator():
    R = np.array([[[0,0],[0,1]]])
    pn = OpenPNM.Network.Template()
    pn.generate(R)
    pn.trim(R>0.5)
    assert len(pn.get_pore_data(prop='coords'))==3
    assert len(pn.get_throat_data(prop='conns'))==2

def test_linear_solver():
    # fix cube dimensions?
    pn = OpenPNM.Network.Cubic()
    pn.generate(add_boundaries=False)

    x,y,z = pn.get_pore_data(prop='coords').T
    left = x==x.min()
    right = x==x.max()
    ics = 2*left + 1*right

    sol = OpenPNM.Shortcuts.solve_linear(pn, ics)
    assert( np.allclose(sol, np.array([2,1.5,1]*9)) )

def test_rectilinear_integrity_after_prune():
    R = np.random.rand(50,40,30)
    # some simple visual pruning, for comparison
    M = np.where(R > R.mean(), R, 0)
    # the convoluted graph way
    pn = OpenPNM.Network.Template()
    pn.generate(R)
    pn.trim(R<=R.mean())
    O = pn.asarray()
    # what it would look like normally
    assert np.allclose(M, O)

if __name__ == '__main__':
    pytest.main([__file__])