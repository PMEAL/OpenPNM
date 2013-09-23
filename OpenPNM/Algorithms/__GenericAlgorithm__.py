#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __GenericAlgorithm__: Base class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np
import scipy.sparse as sprs
from time import clock
import heapq
import itertools

class GenericAlgorithm(OpenPNM.Base.OpenPNMbase):
    r"""
    GenericAlgorithm - Base class to execute algorithms
    
    Parameters
    ----------
    net : Descendent of OpenPNM.Network.GenericNetwork
        A valid network for this algorithm
    
    
    
            
    Examples
    --------
    
    To reserve space for a network with the default number of pores
    and throats execute
    
    >>> import OpenPNM as PNM
    >>> net=PNM.Algorithms.GenericAlgorithm()
    
    
    .. note:: In development
    """
    
    def __init__(self,net=OpenPNM.Network.GenericNetwork,**kwords):
        r"""
        Initialize
        """
        super(GenericAlgorithm,self).__init__(**kwords)
        self.indent = ""
        self._logger.debug("Construct class")
        self._net = net
        
    def run(self):
        r"""
        Main run command for the algorithm
        """
        self._logger.info(self.indent+"Execute run(): Basic version")
        indent=self.indent
        self.indent=self.indent + '  '
        self._do_outer_iteration_stage()
        self.indent=indent
        print self.indent, "-"*39
    
    def _do_one_outer_iteration(self):
        r"""
        One iteration of an outer iteration loop for an algorithm 
        (e.g. time or parametric study)
        """
        self._logger.info(self.indent+"One Outer Iteration: Basic version")
        indent=self.indent
        self.indent=self.indent + '  '
        self._do_inner_iteration_stage()
        self.indent=indent
        print self.indent, "-"*39
        
    def _do_outer_iteration_stage(self):
        r"""
        Executes the outer iteration stage
        """
        self._logger.info(self.indent+"Outer Iteration Stage: Basic version")
        indent=self.indent
        self.indent=self.indent + '  '
        self._do_one_outer_iteration()
        self.indent=indent
        print self.indent, "-"*39
    def _do_one_inner_iteration(self):
        r"""
        Executes one inner iteration
        """
        self._logger.warning(self.indent+"One Inner Iteration: Implement me")
    
    def _do_inner_iteration_stage(self):
        r"""
        Executes the inner iteration stage
        """
        self._logger.info(self.indent+"Inner Iteration Stage: Basic version")
        indent=self.indent
        self.indent=self.indent + '  '
        self._do_one_inner_iteration()
        self.indent=indent
        print self.indent, "-"*39
        

if __name__ =="__main__":
    test = GenericAlgorithm(loggername="TestGenericAlg")
    test.run()