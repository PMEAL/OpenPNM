# -*- coding: utf-8 -*-

import sys, os
import OpenPNM
from OpenPNM.Phases.__GenericPhase__ import GenericPhase
from OpenPNM.Phases import models as fm

class TestPhase(GenericPhase):
    r'''
    Creates Phase object with a default name 'testphase' and preset values for an air-like 
    
    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this phase object will be attached.  
        
    Notes
    -----
    This explicit association is necessary so the Phase object can initialize
    data arrays of the correct size to store network data.
    
    Examples
    --------
    >>> pn = OpenPNM.Network.TestNet()
    >>> water = OpenPNM.Phases.Water(network=pn)
    '''
    def __init__(self,name=None,**kwargs):
        super(TestPhase,self).__init__(name=name,**kwargs)
        self._logger.debug("Construct class")
        self._generate()
        
    def _generate(self):
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0
        self['pore.diffusivity'] = 5.4785e-06
        self['pore.molar_density'] = 40.89
        self['pore.viscosity'] = 1.8443e-05
        self['pore.molecular_weight'] = 28.97
        self['pore.critical_temperature'] = 132.5
        self['pore.critical_volume'] = 0.0883


if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    water = OpenPNM.Phases.Water(network=pn)
