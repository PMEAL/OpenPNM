#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __GenericTopology__: Base class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Topology.__init__.py' file.

"""

import OpenPNM
import scipy as sp

class GenericTopology(OpenPNM.Base.Utilities):
    r"""
    GenericTopology - Base class to construct pore networks

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
        super(GenericTopology,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
        
    def generate(self, network, **params):
        r"""
        Generate the network
        """
#        self._logger.debug("self.generate()")
        self._net = network
        self.generate_setup(**params)
        self.generate_pores()
        self.generate_throats()
#        self.add_boundaries()
        self._logger.debug("end of self.generate()")
        del self._net
        return self

    def generate_setup(self,**params):
        r"""
        Perform applicable preliminary checks and calculations required for generation

        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology.
        """
#        self._logger.error("generation_setup: not implemented")

    def generate_pores(self):
        r"""
        Generate the pores (numbering, types and coordinates)

        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology.
        """
#        self._logger.error("generate_pores: not implemented")

    def generate_throats(self):
        r"""
        Generate the throats (numbering and types)

        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology.
        """
#        self._logger.error("generate_throats: not implemented")

    def add_boundaries(self):
        r"""
        Add boundary pores around network (numbering and types), and create necessary throats

        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology.
        """
#        self._logger.error("add_boundaries: not implemented")

if __name__ == '__main__':
    test=GenericTopology(loggername="TestGenerator")

