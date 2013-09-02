#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __GenericPhysics__: Base class to define pore scale physics
==================================================================

.. warning:: The classes of this module should be loaded through the 'PHYS.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np

class GenericPhysics(OpenPNM.BAS.OpenPNMbase):
    r"""


    """
    
    def __init__(self,net=OpenPNM.NET.GenericNetwork,**kwords):
        r"""
        Initialize
        """
        super(GenericPhysics,self).__init__(**kwords)
        self.indent = ""
        self._logger.debug("Construct class")
        self._net = net
        
    def Washburn(self):
        r'''
        this uses the Washburn equation to relate pore size to entry pressure
        '''
        self._net.throat_properties['Pc_entry'] = -4*0.072*np.cos(np.radians(105))/self._net.throat_properties['diameter']
        

if __name__ =="__main__":
    test = GenericPhysics(loggername="TestGenericPhys")
    test.run()