#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jeff Gostick (jeff@gostick.ca)
# License: TBD
# Copyright (c) 2013


"""
module __Tortuosity__: Network Tortuosity Algorithm
========================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms.__init__.py' file.

"""

import scipy as sp
import numpy as np
import scipy.sparse as sprs
import matplotlib.pyplot as plt

from .__GenericAlgorithm__ import GenericAlgorithm

class Tortuosity(GenericAlgorithm):
    r"""
    Determines the tortuosity of the network using a shortest path search algorithm.

    Parameters
    ----------

    loglevel : integer, optional
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    loggername : string, optional
        Set the name of the logger to be output on the console. Defaults to class name.

    Note
    ----
    n/a
    
    """

    def __init__(self, **kwargs):
        r"""

        """
        super(Tortuosity,self).__init__(**kwargs)
        self._logger.debug("Create Tortuosity Object")
        
    def setup(self,
              inlets,
              outlets,
              **params):
        r'''
        '''
        # Parse args
        self._in = inlets
        self._out = outlets

    def run(self):
        r'''
        '''
        self._net.create_adjacency_matrix()
        graph = pn._adjacency_matrix['csr']['conns']
        temp = spgr.shortest_path(csgraph = graph, method='D', directed = False)
        
        
        




if __name__ == '__main__':
    print('no tests yet')
