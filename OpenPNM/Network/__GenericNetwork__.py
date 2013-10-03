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
    
    Attributes
    ----------
    
    self.pore_properties : dictionary (string, ndarray)
        dictionary containing all pore properties.
        The following properties are created automatically
            - numbering    (int)
    self.throat_properties : dictionary (string, ndarray)
        dictionary containing all throat properties.
        The following properties are created automatically
            - numbering     (int)
            - connections   (int,int)   random integers
    self._num_pores : int
        Number of pores
    self._num_throats : int
        Number of throats
    self._needs_update : bool
        flag if the things need to be updated.
    
    

    
    I have included the class attributes on the same level as
    the devs. This had a weird side effect:
      * The first class instance works fine
      * the second takes the same size as the first.
    I have corrected this by initializing the everyting in the constructor
    again. This seems to have solved the problem, but I am not sure why.
    """
    

    def __init__(self,**kwords):
        r'''
        This is the abstract constructor of the basic network class.  
        
        '''
        
        super(GenericNetwork,self).__init__(**kwords)
        self._logger.debug("Method: Constructor")

        self.pore_properties = {}
        self.throat_properties = {}

        self._logger.info("- Creating default pore properties")
        
        self._logger.info("Constructor completed")

    def create_adjacency_matrix(self,tprop='none',sprsfmt='all',dropzeros=True,diag=False,sym=True):
        r"""

        Generates adjacency matricies in various sparse storage formats

        Parameters
        ----------
        tprop : String
            The throat property to enter into the i,j locations
        sprsfmt : String
            The sparse storage format to use
        dropzeros : Boolean
            Remove 0 elements from tprop, instead of creating 0-weighted link
        diag : Boolean
            blah
        sym : Boolean
            Makes the matrix symmetric about the diagonal

        Returns
        -------
        adj_mat : sparse_matrix, optional
            Returns adjacency matrix in specified format for private use
        
        Notes
        -----
        This can return the specified sparse matrix, but will always write the generated matrix to the network object
        
        Examples
        --------
        >>> print 'nothing yet'
        """
        self._logger.debug('create_adjacency_matrix: Start of method')
        Np   = self.get_num_pores()
        Nt   = self.get_num_throats()
        
        if tprop == 'none':
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
        
        if diag:
            print 'not implimenented yet'
        
        if sym:
            row  = sp.append(row,conn[:,1])
            col  = sp.append(col,conn[:,0])
            data = sp.append(data,data)
        
        self._adjmatrix = sprs.coo_matrix((data,(row,col)),(Np,Np))
        if sprsfmt == 'coo' or sprsfmt == 'all':
            self._adjmatrix._coo = self._adjmatrix
            self._adjmatrix_coo = self._adjmatrix
        if sprsfmt == 'csr' or sprsfmt == 'all':
            self._adjmatrix._csr = self._adjmatrix.tocsr()
            self._adjmatrix_csr = self._adjmatrix.tocsr()
        if sprsfmt == 'lil' or sprsfmt == 'all':
            self._adjmatrix._lil = self._adjmatrix.tolil()
            self._adjmatrix_lil = self._adjmatrix.tolil()
        
    def create_incidence_matrix(self,tprop='none',sprsfmt='all',dropzeros=True):
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
        
        if tprop == 'none':
            dataset = np.ones(Nt)
        else:
            dataset = self.throat_properties[tprop]
        
        if dropzeros:
            ind = dataset>0
        else:
            ind = np.ones_like(dataset,dtype=bool)
        conn = self.throat_properties["connections"][ind]
        row  = conn[:,0]
        row = np.append(row,conn[:,1])
        col = self.throat_properties['numbering'][ind]
        col = np.append(col,col)
        data = np.append(dataset[ind],dataset[ind])
        self._incmatrix = sprs.coo.coo_matrix((data,(row,col)),(Np,Nt))
        if sprsfmt == 'coo' or sprsfmt == 'all':
            self._incmatrix_coo = self._incmatrix
            self._incmatrix._coo = self._incmatrix
        if sprsfmt == 'csr' or sprsfmt == 'all':
            self._incmatrix_csr = self._incmatrix.tocsr()
            self._incmatrix._csr = self._incmatrix.tocsr()
        if sprsfmt == 'lil' or sprsfmt == 'all':
            self._incmatrix_lil = self._incmatrix.tolil()
            self._incmatrix._lil = self._incmatrix.tolil()

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

    def get_connected_pores(self,Tnums=[],flatten=True):
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
        Ps = self.throat_properties["connections"][Tnums]
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
            neighborPs = self._adjmatrix._lil.rows[[Pnums]]
        except:
            self.create_adjacency_matrix() 
            neighborPs = self._adjmatrix_lil.rows[[Pnums]]
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
            neighborTs = self._incmatrix._lil.rows[[Pnums]]
        except:
            self.create_incidence_matrix(sprsfmt='lil')
            neighborTs = self._incmatrix._lil.rows[[Pnums]]
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

    def get_neighbor_pores_props(self,Pnum,flatten=True):
        r"""
        Nothing yet, but this will return the specified property rather than
        just the ID numbers
        
        TODO: Impliment
        """
        
    def get_neighbor_throat_props(self,Pnums,Ttype=[0,1,2,3,4,5,6],flatten=True):
        r"""
        Nothing yet, but this will return the specified property rather than
        just the ID numbers
        
        TODO: Impliment
        """

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
        print "="*50
        print "= Overview of network properties"
        print "-"*50
        print "Basic properties of the network"
        print " - Number of pores:   ", self.get_num_pores()
        print " - Number of throats: ", self.get_num_throats()

        print "Pore properties:"
        for key in self.pore_properties:
            print "\t", key,"\t", self.pore_properties[key].dtype, "\t", self.pore_properties[key].shape

        print "Throat properties:"
        for key in self.throat_properties:
            print "\t", key,"\t", self.throat_properties[key].dtype, "\t", self.throat_properties[key].shape

        print "-"*50
        
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
