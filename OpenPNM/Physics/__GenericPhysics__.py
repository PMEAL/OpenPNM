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

class GenericPhysics(OpenPNM.Utilities.OpenPNMbase):
    r"""


    """
    
    def __init__(self,net=OpenPNM.Network.GenericNetwork,**kwords):
        r"""
        Initialize
        """
        super(GenericPhysics,self).__init__(**kwords)
        self.indent = ""
        self._logger.debug("Construct class")
        self._net = net
        
if __name__ =="__main__":
    test = GenericPhysics(loggername="TestGenericPhys")
    test.run()