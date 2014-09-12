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
        self['pore.molecular_weight'] = 28.96               # kg/kmole
        self['pore.critical_pressure'] = 3.786E6            # Pascal
        self['pore.critical_temperature'] = 132.5           # Kelvin
        self['pore.critical_volume'] = 0.002917             # kg/m3
        self['pore.contact_angle'] = 110.0                  # Degree 
        self.add_model(propname='pore.density',
                       model=fm.density.IdealGas)           # kg/m3
        self.add_model(propname='pore.molar_density',
                       model=fm.molar_density.MolarDensity) # kmole/m3
        self['pore.diffusivity'] = 5.4785E-6                # m2/s
        self['pore.surface_tension'] = 0                    # N/m
        self.add_model(propname='pore.thermal_conductivity',
                       model=fm.thermal_conductivity.AirConductivity) # W/m.K 
        self.add_model(propname='pore.viscosity',
                       model=fm.viscosity.AirViscosity)     # kg/m.s
                       
if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    air = OpenPNM.Phases.Air(network=pn)
    