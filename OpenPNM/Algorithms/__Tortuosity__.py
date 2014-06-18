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
import scipy.sparse.csgraph as spgr
import OpenPNM

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
        
    def estimate_time(self):
        r'''
        '''
        pn_temp = OpenPNM.Network.TestNet()
        graph = pn_temp.create_adjacency_matrix(sprsfmt='csr')
        self._net.tic()
        path = spgr.shortest_path(csgraph = graph, method='D', directed = False)
        t = self._net.toc(quiet=True)
        N = 125
        k = 6
        O = t/(N*(N*k + N*sp.log(N)))
        N = self._net.num_pores()
        k = sp.median(self._net.num_neighbors(pores=self._net.pores()))
        t_est = O*(N*(N*k + N*sp.log(N)))
        print("Based on the network size and PC performance, this algorithm will require: ",t_est,' seconds')
        return 

    def run(self,fluid=None):
        r'''
        '''
        self._logger.warning("This algorithm can take some time. \
        Use the estimate_time method to check time requirements for present network.")
        graph = self._net.create_adjacency_matrix(data=self._net['throat.length'],sprsfmt='csr')
        
        if fluid != None:
            self._fluid = fluid
            if 'throat.occupancy' in self._fluid.props():
                temp = self._net['throat.length']*(self._fluid['throat.occupancy']==1)
                graph = self._net.create_adjacency_matrix(data=temp,sprsfmt='csr',prop='temp')
        
        self._net.tic()
        path = spgr.shortest_path(csgraph = graph, method='D', directed = False)
        self._net.toc()
        
        Px = sp.array(self._net['pore.coords'][:,0],ndmin=2)
        Py = sp.array(self._net['pore.coords'][:,1],ndmin=2)
        Pz = sp.array(self._net['pore.coords'][:,2],ndmin=2)
        
        Cx = sp.square(Px.T - Px)
        Cy = sp.square(Py.T - Py)
        Cz = sp.square(Pz.T - Pz)
        Ds = sp.sqrt(Cx + Cy + Cz)
        
        temp = path/Ds
        
        temp[sp.isnan(temp)] = 0
        temp[sp.isinf(temp)] = 0
        
        return temp
        
        
        




if __name__ == '__main__':
    print('no tests yet')
