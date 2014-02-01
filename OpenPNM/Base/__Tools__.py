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
        
    #--------------------------------------------------------------------------
    '''Setter and Getter Methods'''
    #--------------------------------------------------------------------------
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
        try: phase = self.find_object_by_name(phase) #allow passing of fluid name by string
        except: pass #Accept object
        if phase and not subdomain: getattr(phase,'_'+element+'_data')[prop] = sp.array(data,ndmin=1) #Set fluid property
        elif subdomain and not phase: #Set geometry property
            ind = getattr(self,'get_'+element+'_info')(subdomain)
            try: getattr(self,'_'+element+'_data')[prop] #Test existance of prop
            except: getattr(self,'_'+element+'_data')[prop] = sp.zeros((getattr(self,'get_num_'+element+'s')(),))*sp.nan
            if sp.shape(ind) == sp.shape(data): getattr(self,'_'+element+'_data')[prop][ind] = data
            else: print('data is the wrong size!')
        elif phase and subdomain: #Set pore/throat scale physics property
            ind = getattr(self,'get_'+element+'_info')(subdomain)
            try: getattr(phase,'_'+element+'_data')[prop]
            except: getattr(phase,'_'+element+'_data')[prop] = sp.zeros((getattr(self,'get_num_'+element+'s')(),))
            getattr(phase,'_'+element+'_data')[prop][ind] = sp.array(data,ndmin=1)
        elif not (phase or subdomain): getattr(self,'_'+element+'_data')[prop] = sp.array(data,ndmin=1) #Set topology property

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
        try: subdomain = subdomain.name #allow passing of geometry objects
        except: pass #Otherwise, accept string
        try: phase = self.find_object_by_name(phase) #allow passing of fluid name by string
        except: pass #Accept object
        if phase and not subdomain:
            try: return getattr(phase,'_'+element+'_data')[prop] #Get fluid prop
            except: self._logger.error(phase.name+' does not have the requested '+element+' property: '+prop)           
        elif subdomain and not phase: #Get geometry property
            ind = getattr(self,'get_'+element+'_info')(subdomain)
            try: return getattr(self,'_'+element+'_data')[prop][ind]
            except: self._logger.error(subdomain+' does not have the requested '+element+' property: '+prop)            
        elif phase and subdomain: #Get physics property
            ind = getattr(self,'get_'+element+'_info')(subdomain)
            try: return getattr(phase,'_'+element+'_data')[prop][ind] 
            except: self._logger.error(phase.name+'/'+subdomain+' does not have the requested '+element+' property: '+prop) 
        elif not (phase or subdomain): #Get topology property  
            try: return getattr(self,'_'+element+'_data')[prop]
            except: self._logger.error('Network does not have the requested '+element+' property: '+prop)      
 
    def set_pore_data(self,subdomain='',phase='',prop='',data=''):
        r"""
        Deprecated: See set_data
        """
        self.set_data(element='pore',subdomain=subdomain,phase=phase,prop=prop,data=data)
        
    def get_pore_data(self,subdomain='',phase='',prop=''):
        r"""
        Deprecated: See get_data
        """
        return self.get_data(element='pore',subdomain=subdomain,phase=phase,prop=prop)

    def set_throat_data(self,subdomain='',phase='',prop='',data=''):
        r"""
        Deprecated: See set_data
        """
        self.set_data(element='throat',subdomain=subdomain,phase=phase,prop=prop,data=data)         

    def get_throat_data(self,subdomain='',phase='',prop=''):
        r"""
        Deprecated: See get_data
        """
        return self.get_data(element='throat',subdomain=subdomain,phase=phase,prop=prop)     

    def set_info(self,element='',prop='',data='',indices=False):
        r'''
        '''
        if indices:
            try: getattr(self,'_'+element+'_info')[prop]
            except: getattr(self,'_'+element+'_info')[prop] = sp.zeros((getattr(self,'get_num_'+element+'s')(),),dtype=bool)
            getattr(self,'_'+element+'_info')[prop][data] = True
        else:
            getattr(self,'_'+element+'_info')[prop] = sp.array(data,dtype=bool,ndmin=1)

    def get_info(self,element='',prop='',indices=False):
        r'''
        '''
        if indices:
            return sp.where(getattr(self,'_'+element+'_info')[prop]==True)[0]
        else:
            return getattr(self,'_'+element+'_info')[prop]

    def set_pore_info(self,prop='',data='',indices=False):
        r'''
        '''
        self.set_info(element='pore',prop=prop,data=data,indices=indices)

    def get_pore_info(self,prop='',indices=False):
        r'''
        '''
        return self.get_info(element='pore',prop=prop,indices=indices)
        
    def set_throat_info(self,prop='',data='',indices=False):
        r'''
        '''
        self.set_info(element='throat',prop=prop,data=data,indices=indices)
        
    def get_throat_info(self,prop='',indices=False):
        r'''
        '''
        return self.get_info(element='throat',prop=prop,indices=indices)
        
    #--------------------------------------------------------------------------
    '''Object query methods'''
    #--------------------------------------------------------------------------
    def get_num_pores(self,subdomain=['all']):
        r"""
        Returns the number of pores of the specified subdomain

        Parameters
        ----------
        subdomain : list of strings
            Types of pores to count, defaults to 'all'
            
        Returns
        -------
        Np : int
            Returns the number of pores of the specified type

        """
        #convert string to list, if necessary
        if type(subdomain) == str: subdomain = [subdomain]
        Np = sp.shape(self.get_pore_data(prop='numbering'))[0]
        #Count number of pores of specified type
        if subdomain == ['all']: #return all pores
            return Np
        else:
            temp = sp.zeros((Np,),dtype=bool)
            for item in subdomain: #iterate over subdomain list and accumulate Trues
                temp = temp + self.get_info(prop=item,element='pore')
            return sp.sum(temp) #return sum of Trues

    def get_num_throats(self,subdomain=['all']):
        r"""
        Return the number of throats of the specified subdomain

        Parameters
        ----------

        Returns
        -------
        Nt : int

        """
        #convert string to list, if necessary
        if type(subdomain) == str: subdomain = [subdomain]
        Nt = sp.shape(self.get_throat_data(prop='numbering'))[0]
        #Count number of pores of specified type
        if subdomain == ['all']: #return all pores
            return Nt
        else:
            temp = sp.zeros((Nt,),dtype=bool)
            for item in subdomain: #iterate over subdomain list and accumulate Trues
                temp = temp + self.get_info(prop=item,element='throat')
            return sp.sum(temp) #return sum of Trues

    def get_pore_indices(self,subdomain=['all'],indices=True,mode='union'):
        r'''
        '''
        if type(subdomain) == str: subdomain = [subdomain] #convert string to list, if necessary
        if subdomain == ['all']: #Return full index; easier than get_data(prop='nums')
            if indices:
                ind = sp.r_[0:self.get_num_pores()]
            else:
                ind = sp.ones((self.get_num_pores(),),dtype=bool)
        else:
            if mode == 'union':
                union = sp.zeros((self.get_num_pores(),),dtype=bool)
                for item in subdomain: #iterate over subdomain list and collect all indices
                    union = union + self.get_info(element='pore',prop=item)
                ind = union
            elif mode == 'intersection':
                intersect = sp.ones((self.get_num_pores(),),dtype=bool)
                for item in subdomain: #iterate over subdomain list and collect all indices
                    intersect = intersect*self.get_info(element='pore',prop=item)
                ind = intersect
            if indices: ind = sp.where(ind==True)[0]
        return ind

    def get_throat_indices(self,subdomain=['all'],indices=True,mode='union'):
        r'''
        '''
        if type(subdomain) == str: subdomain = [subdomain] #convert string to list, if necessary
        if subdomain == ['all']: #Return full index; easier than get_data(prop='nums')
            if indices:
                ind = sp.r_[0:self.get_num_throats()]
            else:
                ind = sp.ones((self.get_num_throats(),),dtype=bool)
        else:
            if mode == 'union':
                union = sp.zeros((self.get_num_throats(),),dtype=bool)
                for item in subdomain: #iterate over subdomain list and collect all indices
                    union = union + self.get_info(element='throat',prop=item)
                ind = union
            elif mode == 'intersection':
                intersect = sp.ones((self.get_num_throats(),),dtype=bool)
                for item in subdomain: #iterate over subdomain list and collect all indices
                    intersect = intersect*self.get_info(element='throat',prop=item)
                ind = intersect
            if indices: ind = sp.where(ind==True)[0]
        return ind
        
    def find_object_by_name(self,name):
        for item in self._instances:
            if item.name == name:
                obj = item
        return obj

if __name__ == '__main__':
    test1=GenericNetwork(loggername='Test1')
    test2=GenericNetwork(loglevel=20,loggername='Test2')
