# -*- coding: utf-8 -*-
"""
===============================================================================
module __Physics__: Base class for mananging pore-scale Physics properties
===============================================================================

"""
from OpenPNM.Base import logging
logger = logging.getLogger()
from OpenPNM.Network import GenericNetwork
from OpenPNM.Phases import GenericPhase
import OpenPNM.Physics.models
import scipy as sp

class GenericPhysics(OpenPNM.Base.Core):
    r"""
    Generic class to generate Physics objects

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this Physics should be attached

    phase : OpenPNM Phase object
        The Phase object to which this Physics applies

    pores and/or throats : array_like
        The list of pores and throats where this physics applies. If either are
        left blank this will apply the physics nowhere.  The locations can be
        change after instantiation using ``set_locations()``.

    name : str, optional
        A unique string name to identify the Physics object, typically same as
        instance name but can be anything.  If left blank, and name will be
        generated that include the class name and a random string.

    """

    def __init__(self,network=None,phase=None,pores=[],throats=[],**kwargs):
        super(GenericPhysics,self).__init__(**kwargs)
        logger.name = self.name

        #Associate with Network
        if network is None:
            self._net = GenericNetwork()
        else:
            self._net = network  # Attach network to self
            self._net._physics.append(self)  # Register self with network

        #Associate with Phase
        if phase is None:
            self._phases.append(GenericPhase())
        else:
            phase._physics.append(self)  # Register self with phase
            self._phases.append(phase)  # Register phase with self

        #Initialize a label dictionary in the associated fluid
        self._phases[0]['pore.'+self.name] = False
        self._phases[0]['throat.'+self.name] = False
        self._net['pore.'+self.name] = False
        self._net['throat.'+self.name] = False
        self.set_locations(pores=pores,throats=throats)

    def __getitem__(self,key):
        if key.split('.')[-1] == self.name:
            element = key.split('.')[0]
            return self[element+'.all']
        else:
            return super(GenericPhysics,self).__getitem__(key)

    def set_locations(self,pores=[],throats=[]):
        r'''
        This method can be used to set the pore and throats locations of an
        *empty* object.  Once locations have been set they can not be changed.

        Parameters
        ----------
        pores and throats : array_like
            The list of pores and/or throats where the object should be applied.

        Notes
        -----
        This method is intended to assist in the process of loading saved
        objects.  Save data can be loaded onto an empty object, then the object
        can be reassociated with a Network manually by setting the pore and
        throat locations on the object.
        '''
        pores = sp.array(pores,ndmin=1)
        throats = sp.array(throats,ndmin=1)
        if len(pores) > 0:
            #Check for existing Geometry in pores
            temp = sp.zeros((self._net.Np,),bool)
            for key in self._phases[0].physics():
                temp += self._phases[0]['pore.'+key]
            overlaps = sp.sum(temp*self._net.tomask(pores=pores))
            if overlaps > 0:
                raise Exception('The given pores overlap with an existing Physics object')
            #Initialize locations
            self['pore.all'] = sp.ones((sp.shape(pores)[0],),dtype=bool)
            #Specify Physics locations in Phase dictionary
            self._phases[0]['pore.'+self.name][pores] = True
            self._net['pore.'+self.name][pores] = True
        if len(throats) > 0:
            #Check for existing Geometry in pores
            temp = sp.zeros((self._net.Nt,),bool)
            for key in self._phases[0].physics():
                temp += self._phases[0]['throat.'+key]
            overlaps = sp.sum(temp*self._net.tomask(throats=throats))
            if overlaps > 0:
                raise Exception('The given throats overlap with an existing Physics object')
            #Initialize locations
            self['throat.all'] = sp.ones((sp.shape(throats)[0],),dtype=bool)
            #Specify Physics locations in Phase dictionary
            self._phases[0]['throat.'+self.name][throats] = True
            self._net['throat.'+self.name][throats] = True

if __name__ == '__main__':
    print('none yet')


