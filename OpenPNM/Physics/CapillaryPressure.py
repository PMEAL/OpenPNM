#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __GenericPhysics__: Base class to define pore scale physics
==================================================================

.. warning:: The classes of this module should be loaded through the 'Physics.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np

def Washburn(network,sigma,theta):
    return -4*sigma*np.cos(np.radians(theta))/network.throat_properties['diameter']
    

