#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __GenericGeometry__: Base class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import scipy.stats as spst
from functools import partial
import numpy as np

class GenericGeometry(OpenPNM.Base.Utilities):
    r"""
    GenericGeometry - Base class to construct pore networks

    This class contains the interface definition for the construction of networks

    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)


    """

    def __init__(self, **kwargs):

        r"""
        Initialize
        """
        super(GenericGeometry,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
       
    def create(self,net,**prms):
        r"""
        Create a geometry object using the supplied parameters
        """
        for key, args in prms.items():
            try:
                function = getattr( getattr(OpenPNM.Geometry, key), args['method'] ) # this gets the method from the file
                preloaded_fn = partial(function, geometry=self, network=net, **args) 
                setattr(self, key, preloaded_fn)
                self._logger.info("Successfully loaded {}.".format(key))
            except AttributeError:
                self._logger.debug("Did not manage to load {}.".format(key))
        self.regenerate()
        return self
        
    def regenerate(self):
        self._logger.info("Refreshing geometry")
        self.pore_seed()
        self.pore_diameter()
        self.pore_volume()
        self.throat_seed()
        self.throat_diameter()
        self.throat_length()
        self.throat_volume()
        

if __name__ == '__main__':
    test=GenericGeometry(loggername="TestGenerator")

