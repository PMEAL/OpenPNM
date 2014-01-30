#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __Tools__: Base class to construct pore network tools
==================================================================

.. warning:: The classes of this module should be loaded through the 'Base.__init__.py' file.

"""

import OpenPNM
import numpy as np
import scipy as sp
import scipy.sparse as sprs
import pprint
import collections
from . import Utilities


class Tools(Utilities):
    r"""
    Tools - Base class to initialize pore network methods

    This class contains the interface definition for the construction of networks

    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    loggername : string (optional)
        Define the logger name to be used on console output. Defaults to class name.

    """

    def __init__(self, **kwargs):

        r"""
        Initialize
        """
        super(Tools,self).__init__(**kwargs)
        self._logger.debug("Construct Network container")
        #Initialize network properties dictionaries
        self._pore_data = {}
        self._pore_info = {}
        self._throat_data = {}
        self._throat_info = {}
        #Initialize fluid, physics, and geometry tracking lists
        self._fluids = []
        self._geometry = []
        self._logger.info("Construction of Network container complete")
        
    #-------------------------------------------------------------------
    '''Generalized pore_data setter and getter methods'''
    #-------------------------------------------------------------------
    
    def set_data(self,element='',subdomain='',phase='',prop='',data=''):
        r"""
        Writes data to fluid or network objects according to input arguments.
        Parameters
        ----------
        prop : string
            Name of property to write
        fluid : string, optional
            Name of fluid to which data is written.  If omitted data is written to network object.
        data : array_like
            Data values to write to object
        """
        try: subdomain = subdomain.name #allow passing of geometry objects
        except: pass #Otherwise, accept string
        if phase and not subdomain: getattr(phase,'_'+element+'_data')[prop] = sp.array(data,ndmin=1) #Set fluid property
        elif subdomain and not phase: #Set geometry property
            ind = getattr(self,'get_'+element+'_indices')(subdomain)
            try: getattr(self,'_'+element+'_data')[prop] #Test existance of prop
            except: getattr(self,'_'+element+'_data')[prop] = sp.zeros((getattr(self,'get_num_'+element+'s')(),))*sp.nan
            if sp.shape(ind) == sp.shape(data): getattr(self,'_'+element+'_data')[prop][ind] = data
            else: print('data is the wrong size!')
        elif phase and subdomain: #Set pore/throat scale physics property
            ind = getattr(self,'get_'+element+'_indices')(subdomain)
            try: getattr(phase,'_'+element+'_data')[prop]
            except: getattr(phase,'_'+element+'_data')[prop] = sp.zeros((getattr(self,'get_num_'+element+'s')(),))
            getattr(phase,'_'+element+'_data')[prop][ind] = sp.array(data,ndmin=1)
        elif not (phase or subdomain): getattr(self,'_'+element+'_data')[prop] = sp.array(data,ndmin=1) #Set topology property


    #-------------------------------------------------------------------
    '''Generalized throat_data setter and getter methods'''
    #-------------------------------------------------------------------


    def set_pore_data(self,subdomain='',phase='',prop='',data=''):
        r"""
        Writes pore data to fluid or network objects according to input arguments.
        Parameters
        ----------
        prop : string
            Name of property to write
        fluid : string, optional
            Name of fluid to which data is written.  If omitted data is written to network object.
        data : array_like
            Data values to write to object
        """
        self.set_data(subdomain,phase,prop,data,element='pore')

    def set_throat_data(self,subdomain='',phase='',prop='',data=''):
        r"""
        Writes throat data to fluid or network objects according to input arguments.
        Parameters
        ----------
        prop : string
            Name of property to write
        fluid : string, optional
            Name of fluid to which data is written.  If omitted data is written to network object.
        data : array_like
            Data values to write to object
        """
        self.set_data(subdomain,phase,prop,data,element='throat')    
    
    def get_data(self,element='',subdomain='',phase='',prop=''):
        r"""
        Retrieves data from fluid or network objects according to input arguments.
        Parameters
        ----------
        prop : string
            Name of property to retrieve.  Requesting property 'all' prints a list of existing properties.
        fluid : string, optional
            Name of fluid from which to retrieve data.  If omitted data is retrieved from network object.
        Returns
        -------
        array_like
            An ndarray containing the requested property data from the specified object
        """            
        if phase and not subdomain:
            try: return getattr(phase,'_'+element+'_data')[prop] #Get fluid prop
            except: self._logger.error(phase+' does not have the requested '+element+' property: '+prop)           
        elif subdomain and not phase: #Get geometry property
            ind = getattr(self,'get_'+element+'_indices')(subdomain)
            try: return getattr(self,'_'+element+'_data')[prop][ind]
            except: self._logger.error(subdomain+' does not have the requested '+element+' property: '+prop)            
        elif phase and subdomain: #Get physics property
            ind = getattr(self,'get_'+element+'_indices')(subdomain)
            try: return getattr(phase,'_'+element+'_data')[prop][ind] 
            except: self._logger.error(phase+'/'+subdomain+' does not have the requested '+element+' property: '+prop) 
        elif not (phase or subdomain): #Get topology property  
            try: return getattr(self,'_'+element+'_data')[prop]
            except: self._logger.error('Network does not have the requested '+element+' property: '+prop)      
    
    def get_pore_data(self,subdomain='',phase='',prop=''):
        r"""
        Retrieves pore data from fluid or network objects according to input arguments.
        Parameters
        ----------
        prop : string
            Name of property to retrieve.  Requesting property 'all' prints a list of existing properties.
        fluid : string, optional
            Name of fluid from which to retrieve data.  If omitted data is retrieved from network object.
        Returns
        -------
        array_like
            An ndarray containing the requested property data from the specified object
        """
        self.get_data(subdomain,phase,prop,element='pore')     


    def get_throat_data(self,subdomain='',phase='',prop=''):
        r"""
        Retrieves throat data from fluid or network objects according to input arguments.
        Parameters
        ----------
        prop : string
            Name of property to retrieve.  Requesting property 'all' prints a list of existing properties.
        fluid : string, optional
            Name of fluid from which to retrieve data.  If omitted data is retrieved from network object.
        Returns
        -------
        array_like
            An ndarray containing the requested property data from the specified object
        """
        self.get_data(subdomain,phase,prop,element='throat')     

    def set_info(self,element='',prop='',data='',indices=False):
        r'''
        '''
        if indices:
            try: getattr(self,'_'+element+'_info')[prop]
            except: getattr(self,'_'+element+'_info')[prop] = sp.zeros((getattr(self,'get_num_'+element+'s')(),),dtype=bool)
            getattr(self,'_'+element+'_info')[prop][data] = True
        else:
            getattr(self,'_'+element+'_info')[prop] = sp.array(data,dtype=bool,ndmin=1)
    
    def set_pore_info(self,prop='',data='',indices=False):
        r'''
        '''
        self.set_info(prop,data,indices,element='pore')
        
    def set_throat_info(self,prop='',data='',indices=False):
        r'''
        '''
        self.set_info(prop,data,indices,element='throat')
        
        
    def get_info(self,element='',prop='',indices=False):
        r'''
        '''
        if indices:
            return sp.where(getattr(self,'_'+element+'_info')[prop]==True)[0]
        else:
            return getattr(self,'_'+element+'_info')[prop]

    def get_pore_info(self,prop='',indices=False):
        r'''
        '''
        self.get_info(prop,indices,element='pore')
        
    def get_throat_info(self,prop='',indices=False):
        r'''
        '''
        self.get_info(prop,indices,element='throat')        


    def get_num_pores(self,subdomain=['all']):
        r"""
        Returns the number of pores of the specified subdomain

        Parameters
        ----------

        Returns
        -------
        Np : int
            Returns the number of pores of the specified type

        """
        #Parse Ptype input argument
        Ptype = self.get_type_definitions(Ptype).number

        #Count number of pores of specified type
        try:
            Np = np.sum(np.in1d(self._pore_data['type'],Ptype))
        except:
            Np = 0
        return Np

    def get_num_throats(self,subdomain=['all']):
        r"""
        Return the number of throats of the specified subdomain

        Parameters
        ----------

        Returns
        -------
        Nt : int

        """
        #Parse Ttype input argument
        Ttype = self.get_type_definitions(Ttype).number

        #Count number of throat of specified type
        try:
            Nt = np.sum(np.in1d(self._throat_data['type'],Ttype))
        except:
            Nt = 0
        return Nt

    def get_pore_indices(self,subdomain=['all']):
        r'''
        '''
        if subdomain == ['all']: #Return full index; easier than get_data(prop='nums')
            ind = sp.r_[0:self.get_num_pores()]
        else:
            if type(subdomain) == str: subdomain = [subdomain] #convert sting to list, if necessary
            ind = []
            for item in subdomain: #iterate over subdomain list and collect all indices
                ind.extend(list(self.get_info(prop=item,indices=True)))
            ind = sp.unique(ind) #remove duplicate indices and sort
        return ind

    def get_throat_indices(self,subdomain=['all']):
        r'''
        '''
        if subdomain == ['all']: #Return full index; easier than get_data(prop='nums')
            ind = sp.r_[0:self.get_num_throats()]
        else:
            if type(subdomain) == str: subdomain = [subdomain] #convert sting to list, if necessary
            ind = []
            for item in subdomain: #iterate over subdomain list and collect all indices
                ind.extend(list(self.get_info(prop=item,indices=True)))
            ind = sp.unique(ind) #remove duplicate indices and sort
        return ind


    def fluids_listing(self):
        r"""
        Prints the names of all fluids objects attached to the network
        """
        for item in self._fluids:
            print(item.name+': ',item)

    def fluids_update(self,name='all'):
        r"""
        Updates ALL properties of specified fluid object attached to the network

        Parameters
        ----------
        name : string (optional)
            The name of fluid to be updated.  An empty string (default) refreshes all fluids.
        """
        for item in self._fluids:
            if (item.name == name) or (name == 'all'):
                item.regenerate()
                self._logger.info('Refreshed '+item.name)

    def geometry_listing(self):
        r"""
        """
        for item in self._geometry:
            print(item.name+': ',item)

    def geometry_update(self,name='all'):
        r"""
        """
        for item in self._geometry:
            if (item.name == name) or (name == 'all'):
                item.regenerate()
                self._logger.info('Refreshed '+item.name)

if __name__ == '__main__':
    test1=GenericNetwork(loggername='Test1')
    test2=GenericNetwork(loglevel=20,loggername='Test2')
