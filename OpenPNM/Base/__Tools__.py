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
    def set_data(self,element='',subdomain='',phase='',prop='',data='',indices=''):
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
        if phase and not subdomain: #Set fluid property
            try: getattr(phase,'_'+element+'_data')[prop]
            except: getattr(phase,'_'+element+'_data')[prop] = sp.zeros((getattr(phase,'get_num_'+element+'s')(),))
            if indices!='': getattr(phase,'_'+element+'_data')[prop][indices] = sp.array(data,ndmin=1)
            else: getattr(phase,'_'+element+'_data')[prop] = sp.array(data,ndmin=1) 
        elif subdomain and not phase: #Set geometry property
            ind = getattr(self,'get_'+element+'_info')(subdomain)
            try: getattr(self,'_'+element+'_data')[prop] #Test existance of prop
            except: getattr(self,'_'+element+'_data')[prop] = sp.zeros((getattr(self,'get_num_'+element+'s')(),))*sp.nan
            if indices!='':
                if (sp.in1d(getattr(self,'get_'+element+'_indices')()[indices],\
                getattr(self,'get_'+element+'_indices')(subdomain))).all():
                    ind_temp = sp.zeros((getattr(self,'get_num_'+element+'s')(),),dtype=bool)
                    ind_temp[indices] = True
                    ind = ind_temp
                else: self._logger.error('Some/all of these indices do not belong to this subdomain!')
            if sp.sum(ind) == sp.shape(data)[0] or sp.shape(data)[0]==1:
                getattr(self,'_'+element+'_data')[prop][ind] = sp.array(data,ndmin=1)
            else: print('data is the wrong size!')
                
        elif phase and subdomain: #Set pore/throat scale physics property
            ind = getattr(self,'get_'+element+'_info')(subdomain)
            try: getattr(phase,'_'+element+'_data')[prop]
            except: getattr(phase,'_'+element+'_data')[prop] = sp.zeros((getattr(phase,'get_num_'+element+'s')(),))
            if indices!='':
                if (sp.in1d(getattr(self,'get_'+element+'_indices')()[indices],\
                getattr(self,'get_'+element+'_indices')(subdomain))).all():
                    ind_temp = sp.zeros((getattr(phase,'get_num_'+element+'s')(),),dtype=bool)
                    ind_temp[indices] = True
                    ind = ind_temp
                else: phase._logger.error('Some/all of these indices do not belong to this subdomain!')
            if sp.sum(ind) == sp.shape(data)[0] or sp.shape(data)[0]==1:
                getattr(phase,'_'+element+'_data')[prop][ind] = sp.array(data,ndmin=1)
            else: print('data is the wrong size!')
        elif not (phase or subdomain):  #Set topology property
            try: getattr(self,'_'+element+'_data')[prop]
            except: getattr(self,'_'+element+'_data')[prop] = sp.zeros_like(data)           
            if indices!='': getattr(self,'_'+element+'_data')[prop][indices] = sp.array(data,ndmin=1)
            else: getattr(self,'_'+element+'_data')[prop] = sp.array(data,ndmin=1)

    def get_data(self,element='',subdomain='',phase='',prop='',indices=''):
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
            try: 
                getattr(phase,'_'+element+'_data')[prop]
                if indices!='':  return getattr(phase,'_'+element+'_data')[prop][indices]
                else: return getattr(phase,'_'+element+'_data')[prop] #Get fluid prop
            except: phase._logger.error(phase.name+' does not have the requested '+element+' property: '+prop)           
        elif subdomain and not phase: #Get geometry property
            ind = getattr(self,'get_'+element+'_info')(subdomain)            
            try: 
                getattr(self,'_'+element+'_data')[prop]                
                if indices!='':
                    if (sp.in1d(getattr(self,'get_'+element+'_indices')()[indices],\
                    getattr(self,'get_'+element+'_indices')(subdomain))).all():
                        ind_temp = sp.zeros((getattr(self,'get_num_'+element+'s')(),),dtype=bool)
                        ind_temp[indices] = True
                        ind = ind_temp
                    else: self._logger.error('Some/all of these indices do not belong to this subdomain!')
                return getattr(self,'_'+element+'_data')[prop][ind]
            except: self._logger.error(subdomain+' does not have the requested '+element+' property: '+prop)            
        elif phase and subdomain: #Get physics property
            ind = getattr(self,'get_'+element+'_info')(subdomain)            
            try: 
                getattr(phase,'_'+element+'_data')[prop]
                if indices!='':
                    if (sp.in1d(getattr(self,'get_'+element+'_indices')()[indices],\
                    getattr(self,'get_'+element+'_indices')(subdomain))).all():
                        ind_temp = sp.zeros((getattr(phase,'get_num_'+element+'s')(),),dtype=bool)
                        ind_temp[indices] = True
                        ind = ind_temp
                    else: phase._logger.error('Some/all of these indices do not belong to this subdomain!')                   
                return getattr(phase,'_'+element+'_data')[prop][ind]
            except: phase._logger.error(phase.name+'/'+subdomain+' does not have the requested '+element+' property: '+prop) 
        elif not (phase or subdomain): #Get topology property  
            try: 
                getattr(self,'_'+element+'_data')[prop]
                if indices!='':  return getattr(self,'_'+element+'_data')[prop][indices]
                else: return getattr(self,'_'+element+'_data')[prop] #Get fluid prop
            except: self._logger.error('Object does not have the requested '+element+' property: '+prop)      
 
    def set_pore_data(self,subdomain='',phase='',prop='',data='',indices=''):
        r"""
        Deprecated: See set_data
        """
        self.set_data(element='pore',subdomain=subdomain,phase=phase,prop=prop,data=data,indices=indices)
        
    def get_pore_data(self,subdomain='',phase='',prop='',indices=''):
        r"""
        Deprecated: See get_data
        """
        return self.get_data(element='pore',subdomain=subdomain,phase=phase,prop=prop,indices=indices)

    def set_throat_data(self,subdomain='',phase='',prop='',data='',indices=''):
        r"""
        Deprecated: See set_data
        """
        self.set_data(element='throat',subdomain=subdomain,phase=phase,prop=prop,data=data,indices=indices)         

    def get_throat_data(self,subdomain='',phase='',prop='',indices=''):
        r"""
        Deprecated: See get_data
        """
        return self.get_data(element='throat',subdomain=subdomain,phase=phase,prop=prop,indices=indices)     

    def set_info(self,element='',prop='',locations='',is_indices=False,mode='merge'):
        r'''
        '''
        if mode=='overwrite':
            getattr(self,'_'+element+'_info')[prop] = sp.zeros((getattr(self,'get_num_'+element+'s')(),),dtype=bool)
        if is_indices:
            try: 
                getattr(self,'_'+element+'_info')[prop]
            except: getattr(self,'_'+element+'_info')[prop] = sp.zeros((getattr(self,'get_num_'+element+'s')(),),dtype=bool)
            getattr(self,'_'+element+'_info')[prop][locations] = True
        else:
            getattr(self,'_'+element+'_info')[prop] = sp.array(locations,dtype=bool,ndmin=1)

    def get_info(self,element='',prop='',return_indices=False):
        r'''
        '''
        if return_indices:
            return sp.where(getattr(self,'_'+element+'_info')[prop]==True)[0]
        else:
            return getattr(self,'_'+element+'_info')[prop]

    def set_pore_info(self,prop='',locations='',is_indices=False):
        r'''
        '''
        self.set_info(element='pore',prop=prop,locations=locations,is_indices=is_indices)

    def get_pore_info(self,prop='',return_indices=False):
        r'''
        '''
        return self.get_info(element='pore',prop=prop,return_indices=return_indices)
        
    def set_throat_info(self,prop='',locations='',is_indices=False):
        r'''
        '''
        self.set_info(element='throat',prop=prop,locations=locations,is_indices=is_indices)
        
    def get_throat_info(self,prop='',return_indices=False):
        r'''
        '''
        return self.get_info(element='throat',prop=prop,return_indices=return_indices)
        
    #--------------------------------------------------------------------------
    '''Object query methods'''
    #--------------------------------------------------------------------------
    def get_num_pores(self,subdomain=['all'],mode='union'):
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
        #Count number of pores of specified type
        if subdomain == ['all']: #return all pores
            return sp.shape(self.get_pore_info(prop='numbering'))[0]
        else:
            temp = self.get_pore_indices(subdomain=subdomain,mode=mode,indices=False)
            return sp.sum(temp) #return sum of Trues
            

    def get_num_throats(self,subdomain=['all'],mode='union'):
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
        #Count number of pores of specified type
        if subdomain == ['all']: #return all pores
            return sp.shape(self.get_throat_info(prop='numbering'))[0]
        else:
            temp = self.get_throat_indices(subdomain=subdomain,mode=mode,indices=False)
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

    def get_result(self,obj,**kwargs):
        
        obj.update(**kwargs)        


if __name__ == '__main__':
    test1=GenericNetwork(loggername='Test1')
    test2=GenericNetwork(loglevel=20,loggername='Test2')
