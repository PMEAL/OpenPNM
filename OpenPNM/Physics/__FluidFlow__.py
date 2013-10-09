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

class FluidFlow(OpenPNM.Utilities.OpenPNMbase):
    r"""
    Methods in this class are used to determine the hydraulic conductivity of pores and throats from their geometric properties.
    """
    def __init__(self,**kwargs):
        super(FluidFlow,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")

    def HagenPoiseuille(self,network,viscosity):
        r"""
        Calculates the hydraulic conductvity of throat assuming cylindrical geometry
        
        Parameters
        ----------
        network : OpenPNM Object
            The network for which the calculations are desired
            
        viscosity : float
            The viscosity of the fluid
        """
        
        return 'nothing yet'

