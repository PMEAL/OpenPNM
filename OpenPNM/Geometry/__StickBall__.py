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
from .__GenericGeometry__ import GenericGeometry

class Stick_and_Ball(GenericGeometry):
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
        super(Stick_and_Ball,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
   
        self.add_method(prop='pore_seed',model='random')
        self.add_method(prop='throat_seed',model='neighbor_min')
        self.add_method(prop='pore_diameter',model='sphere',name='weibull_min',shape=2.5,loc=6e-6,scale=2e-5)
        self.add_method(prop='throat_diameter',model='cylinder',name='weibull_min',shape=2.5,loc=6e-6,scale=2e-5)
        self.add_method(prop='pore_volume',model='sphere')
        self.add_method(prop='throat_length',model='straight')
        self.add_method(prop='throat_volume',model='cylinder')
        

if __name__ == '__main__':
    test=GenericGeometry(loggername="TestGenerator")

