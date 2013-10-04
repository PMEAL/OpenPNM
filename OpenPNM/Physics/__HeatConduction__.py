#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __GenericPhysics__: Base class to define pore scale physics
==================================================================

.. warning:: The module must be loaded by adding an import line to Physics/__init__.py

"""

import OpenPNM
import scipy as sp

class HeatConduction(OpenPNM.Utilities.OpenPNMbase):
    def __init__(self):
        print 'init Multiphase'

    def FourierConductiity(self,network):
        r"""
        ---
        """
        
        return
