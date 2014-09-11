import sys, os
import OpenPNM
from OpenPNM.Phases.__GenericPhase__ import GenericPhase
from OpenPNM.Phases import models as fm

class Air(GenericPhase):
    r"""
    Creates Phase object with a default name 'air' and preset values for air
    
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
    >>> air = OpenPNM.Phases.Air(network=pn)
    """
    def __init__(self,name=None,**kwargs):
        super(Air,self).__init__(name=name,**kwargs)
        self._logger.debug("Construct class")
        self._generate()
        
    def _generate(self):
        self['pore.molecular_weight'] = 28.96       # kg/kmole
        self['pore.critical_pressure'] = 3.786E6    # Pascal
        self['pore.critical_temperature'] = 132.5   # Kelvin
        self['pore.critical_volume'] = 0.002917     # kg/m3
        self['pore.contact_angle'] = 110.0          # Degree 
        self['pore.density'] = 1.185                # kg/m3
        self['pore.diffusivity'] = 5.4785E-6        # m2/s
        self['pore.molar_density'] = 0.0409         # kmole/m3
        self['pore.surface_tension'] = 275.5E-6     # N/m
        self['pore.thermal_conductivity'] = 0.02624 # W/m.K      
        self['pore.viscosity'] = 18.44E-6           # kg/m.s                      
        
if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    air = OpenPNM.Phases.Air(network=pn)
    