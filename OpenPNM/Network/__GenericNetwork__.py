#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __GenericNetwork__: contains OpenPNM topology baseclass
==============================================================


.. warning:: The classes of this module should be loaded through the 'NET.__init__.py' file.

"""

import OpenPNM
import numpy as np
import scipy as sp
import scipy.sparse as sprs
import matplotlib as mpl

class GenericNetwork(OpenPNM.Utilities.OpenPNMbase):
    r"""
    GenericNetwork - Base topology class for pore networks

    This class contains the basic functionality for storaging and querying pore network data.

    Parameters
    ----------
    num_pores : int
        Number of pores
    num_throats : int
        Number of throats
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)

    """

    def __init__(self,**kwords):
        r'''
        This is the abstract constructor of the basic network class.

        '''

        super(GenericNetwork,self).__init__(**kwords)
        self._logger.debug("Method: Constructor")

        #Initializes properties dictionaries
        self.pore_properties = {}
        self.throat_properties = {}

        #Initialize the conditions dictionaries
        self.pore_conditions = {}
        self.throat_conditions = {}

        #Initialize constants dictionary
        self.const = {}

        #This initializes the custom 'self-protecting' dictionary
#        self.pore_properties = {}
#        self.throat_properties = {}

        #Initialize adjacency and incidence matrix dictionaries
        self.adjacency_matrix = {}
        self.incidence_matrix = {}
        self.adjacency_matrix['coo'] = {}
        self.adjacency_matrix['csr'] = {}
        self.adjacency_matrix['lil'] = {}
        self.incidence_matrix['coo'] = {}
        self.incidence_matrix['csr'] = {}
        self.incidence_matrix['lil'] = {}

        self._logger.info("Constructor completed")

    def create_adjacency_matrix(self,tprop='connections',sprsfmt='all',dropzeros=True,sym=True):
        r"""

        Generates adjacency matricies in various sparse storage formats

        Parameters
        ----------
        tprop : String, optional
            The throat property to enter into the i,j locations. If no property is given 'connections' is used.
        sprsfmt : String, optional
            The sparse storage format to use. If none type is given, all are generated (coo, csr & lil)
        dropzeros : Boolean, optional
            Remove 0 elements from tprop, instead of creating 0-weighted links, the default is True.
        sym : Boolean, optional
            Makes the matrix symmetric about the diagonal, the default is true.

        Returns
        -------
        adj_mat : sparse_matrix, optional
            Returns adjacency matrix in specified format for private use.

        Notes
        -----
        This 'can' return the specified sparse matrix, but will always write the generated matrix to the network object

        Examples
        --------
        >>> print 'nothing yet'
        """
        self._logger.debug('create_adjacency_matrix: Start of method')
        Np   = self.get_num_pores()
        Nt   = self.get_num_throats()

        if tprop == 'connections':
            dataset = np.ones(Nt)
        else:
            dataset = self.throat_properties[tprop]

        if dropzeros:
            ind = dataset>0
        else:
            ind = np.ones_like(dataset,dtype=bool)

        conn = self.throat_properties["connections"][ind]
        row  = conn[:,0]
        col  = conn[:,1]
        data = dataset[ind]

        #Append row & col to each other, and data to itself
        if sym:
            row  = sp.append(row,conn[:,1])
            col  = sp.append(col,conn[:,0])
            data = sp.append(data,data)

        temp = sprs.coo_matrix((data,(row,col)),(Np,Np))
        if sprsfmt == 'coo' or sprsfmt == 'all':
            self.adjacency_matrix['coo'][tprop] = temp
        if sprsfmt == 'csr' or sprsfmt == 'all':
            self.adjacency_matrix['csr'][tprop] = temp.tocsr()
        if sprsfmt == 'lil' or sprsfmt == 'all':
            self.adjacency_matrix['lil'][tprop] = temp.tolil()
        if sprsfmt != 'all':
            return self.adjacency_matrix[sprsfmt][tprop]

    def create_incidence_matrix(self,tprop='connections',sprsfmt='all',dropzeros=True):
        r"""

        Creates an incidence matrix filled with specified throat values

        Parameters
        ----------
        tname : The property of interest (i.e. diameter, volume, etc.)

        Returns
        -------
        An incidence matrix (a cousin to the adjacency matrix, useful for finding throats of given a pore)

        Notes
        -----

        Examples
        --------
        >>> print 'nothing yet'
        """
        self._logger.debug('create_incidence_matrix: Start of method')

        Nt = self.get_num_throats()
        Np = self.get_num_pores()

        if tprop == 'connections':
            dataset = np.ones(Nt)
        else:
            dataset = self.throat_properties[tprop]

        if dropzeros:
            ind = dataset>0
        else:
            ind = np.ones_like(dataset,dtype=bool)

        conn = self.throat_properties['connections'][ind]
        row  = conn[:,0]
        row = np.append(row,conn[:,1])
        col = self.throat_properties['numbering'][ind]
        col = np.append(col,col)
        data = np.append(dataset[ind],dataset[ind])

        temp = sprs.coo.coo_matrix((data,(row,col)),(Np,Nt))
        if sprsfmt == 'coo' or sprsfmt == 'all':
            self.incidence_matrix['coo'][tprop] = temp
        if sprsfmt == 'csr' or sprsfmt == 'all':
            self.incidence_matrix['csr'][tprop] = temp.tocsr()
        if sprsfmt == 'lil' or sprsfmt == 'all':
            self.incidence_matrix['lil'][tprop] = temp.tolil()
        if sprsfmt != 'all':
            return self.incidence_matrix[sprsfmt][tprop]

    def get_num_pores(self,Ptype=[0,1,2,3,4,5,6]):
        r"""
        Returns the number of pores of the specified type

        Parameters
        ----------

        Ptype : array_like, optional
            list of desired pore types to count

        Returns
        -------
        Np : int

        """
        try:
            Np = np.sum(np.in1d(self.pore_properties['type'],Ptype))
        except:
            Np = 0
        return Np

    def get_num_throats(self,Ttype=[-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6]):
        r"""
        Return the number of throats of the specified type

        Parameters
        ----------

        Ttype : array_like, optional
            list of desired throat types to count

        Returns
        -------
        Nt : int

        """
        try:
            Nt = np.sum(np.in1d(self.throat_properties['type'],Ttype))
        except:
            Nt = 0
        return Nt

    def get_connected_pores(self,Tnums=[],flatten=False):
        r"""
        Return a list of pores connected to a list of throats

        Parameters
        ----------
        Tnums: array_like
            List of throats ID numbers
        flatten : boolean, optional
            If flatten is True (default) a 1D array of unique pore ID numbers
            is returned. If flatten is False the returned array contains
            arrays of neighboring pores for each input throat, in the order
            they were sent.

        Returns
        -------
        Ps : 1D array (if flatten is True) or ndarray of arrays (if flatten is
            False)

        Examples
        --------
        >>> Tnums = [0,1]
        >>> Ps = pn.get_connected_pores(Tnums)
        >>> Ps
        array([  0,   2, 920])

        >>> Tnums = [0,1]
        >>> Ps = pn.get_connected_pores(Tnums,flatten=False)
        >>> Ps
        array([[  0, 920],
               [  0,   2]])
        """
        Ps = self.throat_properties['connections'][Tnums]
#        Ps = [np.asarray(x) for x in Ps if x]
        if flatten:
            Ps = np.unique(np.hstack(Ps))
        return Ps

    def get_connecting_throat(self,P1,P2):
        r"""
        Return a the throat number connecting two given pores connected

        Parameters
        ----------
        Pnum1 , Pnum2 : int

        Returns
        -------
        Tnum : int
            Returns throat ID number, or empty array if pores are not connected
        """
        return np.intersect1d(self.get_neighbor_throats(P1),self.get_neighbor_throats(P2))

    def get_neighbor_pores(self,Pnums,Ptype=[0,1,2,3,4,5,6],flatten=True):
        r"""
        Returns a list of pores neighboring the given pore(s)

        Parameters
        ----------
        Pnums : array_like
            ID numbers of pores whose neighbors are sought.
        Ptype : array_like
            Type of pores to be returned
        flatten : boolean, optional
            If flatten is True (default) a 1D array of unique pore ID numbers
            is returned with the input pores (Pnum) removed. If flatten is
            False the returned array contains arrays of neighboring pores for
            each input pore, in the order they were sent.

        Returns
        -------
        neighborPs : 1D array (if flatten is True) or ndarray of ndarrays (if
            flatten if False)

        Examples
        --------
        >>> Pnums = [0,1]
        >>> Ps = pn.get_neighbor_pores(Pnums)
        >>> Ps
        array([  2,   3, 920, 921])

        >>> Pnums = [0,1]
        >>> Ps = pn.get_neighbor_pores(Pnums,flatten=False)
        >>> Ps
        array([[  1,   2, 920],
               [  0,   3, 921]])
        """
        try:
            neighborPs = self.adjacency_matrix['lil']['connections'].rows[[Pnums]]
        except:
            self.create_adjacency_matrix()
            neighborPs = self.adjacency_matrix['lil']['connections'].rows[[Pnums]]
        if flatten:
            #All the empty lists must be removed to maintain data type after hstack (numpy bug?)
            neighborPs = [sp.asarray(x) for x in neighborPs if x]
            neighborPs = sp.hstack(neighborPs)
            #Remove references to input pores and duplicates
            neighborPs = sp.unique(neighborPs[~np.in1d(neighborPs,Pnums)])
            #Remove pores of the wrong type
            neighborPs = neighborPs[sp.in1d(self.pore_properties['type'][neighborPs],Ptype)]
        else:
            for i in range(0,sp.shape(Pnums)[0]):
                ans = sp.array(sp.where(sp.in1d(self.pore_properties['type'][neighborPs[i]],Ptype)))[0]
                neighborPs[i] = sp.array(neighborPs[i])[ans]
        return np.array(neighborPs)

    def get_neighbor_throats(self,Pnums,Ttype=[0,1,2,3,4,5,6],flatten=True):
        r"""
        Returns a list of throats neighboring the given pore(s)

        Parameters
        ----------
        Pnums : array_like
            ID numbers of pores whose neighbors are sought
        Ttype : array_like
            Type of throats to be returned
        flatten : boolean, optional
            If flatten is True (default) a 1D array of unique throat ID numbers
            is returned. If flatten is False the returned array contains arrays
            of neighboring throat ID numbers for each input pore, in the order
            they were sent.

        Returns
        -------
        neighborTs : 1D array (if flatten is True) or ndarray of arrays (if
            flatten if False)

        Examples
        --------
        >>> Pnums = [0,1]
        >>> Ts = pn.get_neighbor_throats(Pnums)
        >>> Ts
        array([    0,     1,     2, 83895, 83896])

        >>> Pnums = [0,1]
        >>> Ts = pn.get_neighbor_throats(Pnums,flatten=False)
        >>> Ts
        array([[    0,     1,     2],
               [    2, 83895, 83896]])
        """
        try:
            neighborTs = self.incidence_matrix['lil']['connections'].rows[[Pnums]]
        except:
            self.create_incidence_matrix()
            neighborTs = self.incidence_matrix['lil']['connections'].rows[[Pnums]]
        if flatten and neighborTs:
            #All the empty lists must be removed to maintain data type after hstack (numpy bug?)
            neighborTs = [np.asarray(x) for x in neighborTs if x]
            neighborTs = np.unique(np.hstack(neighborTs))
            #Remove throats of the wrong type
            neighborTs = neighborTs[sp.in1d(self.throat_properties['type'][neighborTs],Ttype)]
        else:
            for i in range(0,sp.shape(Pnums)[0]):
                ans = sp.array(sp.where(sp.in1d(self.throat_properties['type'][neighborTs[i]],Ttype)))[0]
                neighborTs[i] = sp.array(neighborTs[i])[ans]
        return np.array(neighborTs)

    def get_num_neighbors(self,Pnums,Ptype=[0,1,2,3,4,5,6]):
        r"""
        Returns an ndarray containing the number of pores for each element in Pnums

        Parameters
        ----------
        Pnums : array_like
            ID numbers of pores whose neighbors are sought
        Ptype : array_like
            Type of throats to be returne

        Returns
        -------
        num_neighbors : 1D array with number of neighbors in each element

        Examples
        --------
        >>> Pnum = [0,1]
        >>> Nn = pn.get_num_neighbors(Pnum)
        >>> Nn
        array([3, 4], dtype=int8)

        >>> Pnum = range(0,pn.get_num_pores())
        >>> Nn = pn.get_num_neighbors(Pnum)
        >>> Nn
        array([3, 4, 4, ..., 4, 4, 3], dtype=int8)
        >>> pn.pore_properties['num_neighbors'] = Nn
        """
        neighborPs = self.get_neighbor_pores(Pnums,Ptype,flatten=False)
        num = sp.zeros(sp.shape(neighborPs),dtype=sp.int8)
        for i in range(0,sp.shape(num)[0]):
            num[i] = sp.size(neighborPs[i])
        return num

    def get_neighbor_pores_props(self,Pprop,Pnums,Ptype=[0,1,2,3,4,5,6],flatten=True):
        r"""
        Return the desired property for the requested pore ID numbers

        """
        neighborPs = self.get_neighbor_pores(Pnums,Ptype,flatten)
        if flatten:
            propPs = self.pore_properties[Pprop][neighborPs]
        else:
            propPs = self.pore_properties[Pprop][neighborPs]
        return propPs

    def get_neighbor_throat_props(self,Pnums,Ttype=[0,1,2,3,4,5,6],flatten=True):
        r"""
        Nothing yet, but this will return the specified property rather than
        just the ID numbers

        TODO: Impliment
        """

    def interpolate_pore_conditions(self,Tcond=None):
        r"""
        Determines a pore property as the average of it's neighboring throats

        Parameters
        ----------
        Tcond : string
            The name of the throat condition to be interpolated

        Notes
        -----
        This uses an unweighted average, without attempting to account for distances or sizes of pores and throats.

        """
        self.pore_conditions[Tcond] = sp.zeros((self.get_num_pores()))
        #Only interpolate conditions for internal pores, type=0
        Pnums = sp.r_[0:self.get_num_pores(Ptype=[0])]
        nTs = self.get_neighbor_throats(Pnums,flatten=False)
        for i in sp.r_[0:sp.shape(nTs)[0]]:
            self.pore_conditions[Tcond][i] = sp.mean(self.throat_conditions[Tcond][nTs[i]])

    def interpolate_throat_conditions(self,Pcond=None):
        r"""
        Determines a throat condition as the average of the conditions it's neighboring pores

        Parameters
        ----------
        Pcond : string
            The name of the throat condition to be interpolated

        Notes
        -----
        This uses an unweighted average, without attempting to account for distances or sizes of pores and throats.

        """
        self.throat_conditions[Pcond] = sp.zeros((self.get_num_throats()))
        #Interpolate values for all throats, including those leading to boundary
        Tnums = sp.r_[0:self.get_num_throats()]
        nPs = self.get_connected_pores(Tnums,flatten=False)
        for i in sp.r_[0:sp.shape(nPs)[0]]:
            self.throat_conditions[Pcond][i] = sp.mean(self.pore_conditions[Pcond][nPs[i]])

    def check_basic(self):
        r"""
        Check the network for general health
        TODO: implement
        """
        self._logger.debug("Method: check for general healts")

    def print_overview(self):
        r"""
        Print some basic properties
        """
        self._logger.debug("Method: print_overview")
        print "="*72
        print "= Overview of network properties"
        print "-"*72
        print "Basic properties of the network"
        print " - Number of pores:   ", self.get_num_pores()
        print " - Number of throats: ", self.get_num_throats()

        print "="*72
        print "Pore Properties"
        print "-"*72
        print 'PROPERTY', "\t", "\t", 'DTYPE', "\t", 'SHAPE', "\t", 'MEMORY [MB]'
        print "-"*72
        for key in self.pore_properties:
            print key, "\t", "\t", self.pore_properties[key].dtype, "\t", self.pore_properties[key].shape, "\t", self.pore_properties[key].nbytes/1e6

        print "="*72
        print "Throat Properties"
        print "-"*72
        print 'PROPERTY', "\t", "\t", 'DTYPE', "\t", 'SHAPE', "\t", 'MEMORY [MB]'
        print "-"*72
        for key in self.throat_properties:
            print key, "\t", "\t", self.throat_properties[key].dtype, "\t", self.throat_properties[key].shape, "\t", self.throat_properties[key].nbytes/1e6

        print "="*72
        print "Adjacency Matrices"
        print "-"*72
        print 'FORMAT', "\t", 'PROPERTIES'
        print "-"*72
        for sprsfmt in self.adjacency_matrix.keys():
            print sprsfmt, ":\t", self.adjacency_matrix[sprsfmt].keys()

        print "="*72
        print "Incidence Matrices"
        print "-"*72
        print 'FORMAT', "\t", 'PROPERTIES'
        print "-"*72
        for sprsfmt in self.incidence_matrix.keys():
            print sprsfmt, ":\t", self.incidence_matrix[sprsfmt].keys()

        print "="*72

    def __str__(self):
        r"""
        Print some basic properties
        """
        self._logger.debug("Method: print_overview")

        str_overview = (
        '==================================================\n'
        'Overview of network properties\n'
        '--------------------------------------------------\n'
        'Basic properties of the network\n'
        '- Number of pores:   {num_pores}\n'
        '- Number of throats: {num_throats}\n'
        ).format(num_pores=self.get_num_pores(),
                   num_throats=self.get_num_throats())

        str_pore = "\nPore properties:"
        for key, value in self.pore_properties.iteritems():
            str_pore += "\n\t{0:20}{1.dtype:20}{1.shape:20}".format(key,value)

        str_throat = "\nThroat properties:"
        for key, value in self.throat_properties.iteritems():
            str_throat += "\n\t{0:20}{1.dtype:20}{1.shape:20}".format(key,value)

        return str_overview+str_pore+str_throat

    def update(self):
        self.create_adjacency_matrix()
        self.create_incidence_matrix()

if __name__ == '__main__':
    test1=GenericNetwork(loggername='Test1')
    test2=GenericNetwork(loglevel=20,loggername='Test2')
