import sys, os
import OpenPNM
from OpenPNM.Phases.__GenericPhase__ import GenericPhase
from OpenPNM.Phases import models as fm

class Water(GenericPhase):
    r'''
    Creates Phase object with a default name 'water' and preset values for water
    
    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this phase object will be attached.  
        
    Notes
    -----
    This explicit association is necessary so the Phase object can initialize
    data arrays of the correct size to store network data.
    The initial properties are all at std conditions of T = 298 K and P = 1 atm.
    
    Examples
    --------
    >>> pn = OpenPNM.Network.TestNet()
    >>> water = OpenPNM.Phases.Water(network=pn)
    '''
    def __init__(self,name=None,**kwargs):
        super(Water,self).__init__(name=name,**kwargs)
        self._logger.debug("Construct class")
        self._generate()
        
    def _generate(self):
        self['pore.molecular_weight'] = 18.02       # kg/kmole
        self['pore.critical_pressure'] = 2.2064E7   # Pascal
        self['pore.critical_temperature'] = 647.1   # Kelvin
        self['pore.critical_volume'] = 0.003106     # kg/m3
        self['pore.contact_angle'] = 110.0          # Degree 
        self['pore.density'] = 997.1                # kg/m3
        self['pore.diffusivity'] = 1e-9             # m2/s
        self['pore.molar_density'] = 55.35          # kmole/m3
        self['pore.surface_tension'] = 0.07199      # N/m
        self['pore.thermal_conductivity'] = 0.5945  # W/m.K      
        self['pore.vapor_pressure'] = 3141          # Pascal
        self['pore.viscosity'] = 8.936E-4           # kg/m.s

if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    water = OpenPNM.Phases.Water(network=pn)
