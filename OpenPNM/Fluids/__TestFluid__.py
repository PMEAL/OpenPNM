# -*- coding: utf-8 -*-

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM
from OpenPNM.Fluids.__GenericFluid__ import GenericFluid
from OpenPNM.Fluids import models as fm

class TestFluid(GenericFluid):
    r'''
    Creates Fluid object with a default name 'testfluid' and preset values for an air-like 
    
    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this fluid object will be attached.  
        
    Notes
    -----
    This explicit association is necessary so the Fluid object can initialize
    data arrays of the correct size to store network data.
    
    Examples
    --------
    >>> pn = OpenPNM.Network.TestNet()
    >>> water = OpenPNM.Fluids.Water(network=pn)
    '''
    def __init__(self,name=None,**kwargs):
        super(TestFluid,self).__init__(name=name,**kwargs)
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
    water = OpenPNM.Fluids.Water(network=pn)
