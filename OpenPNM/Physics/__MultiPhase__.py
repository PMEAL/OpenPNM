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

class MultiPhase(OpenPNM.Utilities.OpenPNMbase):
    r"""
    Methods in this class are used to track and calculate any pore and throat properties that are dependent on the phase and extent of filling.
    """
    
    def __init__(self,**kwargs):
        super(MultiPhase,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")

    def Conduit_Filled_State_Calculator(self,network):
        r"""
        ---
        """
        
        return

    def Apply_Phase_State_to_Conduit_Conductivity(self):
        r"""
        ---
        """
        
        return 

    def Late_Pore_Filling_Tracking(self,network):
        r"""
        ---
        """

        return