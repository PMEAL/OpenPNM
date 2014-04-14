#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Tests for OpenPNM """

from __future__ import print_function, division, absolute_import
import pytest

import OpenPNM
import numpy as np

def test_template_generator():
    R = np.array([[[0,0],[0,1]]])
    pn = OpenPNM.Network.Template(name='template')
    pn.generate(R, threshold=0.5)
    assert len(pn.get_pore_data(prop='coords'))==3
    assert len(pn.get_throat_data(prop='connections'))==2

if __name__ == '__main__':
    pytest.main([__file__])