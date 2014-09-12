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
        self['pore.molecular_weight'] = 18.02               # kg/kmole
        self['pore.critical_pressure'] = 2.2064E7           # Pascal
        self['pore.critical_temperature'] = 647.1           # Kelvin
        self['pore.critical_volume'] = 0.003106             # kg/m3
        self['pore.contact_angle'] = 110.0                  # Degree 
        self.add_model(propname='pore.density',
                       model=fm.density.WaterDensity)       # kg/m3
        self.add_model(propname='pore.molar_density',
                       model=fm.molar_density.MolarDensity) # kmole/m3
        self['pore.diffusivity'] = 1e-9                     # m2/s
        self.add_model(propname='pore.surface_tension',
                       model=fm.surface_tension.WaterSurfaceTension) # N/m
        self.add_model(propname='pore.thermal_conductivity',
                       model=fm.thermal_conductivity.WaterConductivity) # W/m.K
        self['pore.vapor_pressure'] = 3141                  # Pascal
        self.add_model(propname='pore.viscosity',
                       model=fm.viscosity.WaterViscosity)   # kg/m.s

if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    water = OpenPNM.Phases.Water(network=pn)
