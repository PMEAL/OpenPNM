import sys, os
import OpenPNM
from OpenPNM.Phases.__GenericPhase__ import GenericPhase
from OpenPNM.Phases import models as fm

class Mercury(GenericPhase):
    r"""
    Creates Phase object with a default name 'Hg' and preset values for mercury
    
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
    >>> hg = OpenPNM.Phases.Mercury(network=pn)
    """
    def __init__(self,**kwargs):
        super(Mercury,self).__init__(name='Hg',**kwargs)
        self._logger.debug("Construct class")
        self._generate()
        
    def _generate(self):
        self['pore.molecular_weight'] = 200.6       # kg/kmole
        self['pore.critical_pressure'] = 1.662E8    # Pascal
        self['pore.critical_temperature'] = 1733    # Kelvin
        self['pore.critical_volume'] = 0.000189     # kg/m3
        self['pore.contact_angle'] = 140            # Degree 
        self['pore.density'] = 13534                # kg/m3
        self['pore.diffusivity'] = 1e-9             # m2/s
        self['pore.molar_density'] = 67.467         # kmole/m3
        self['pore.surface_tension'] = 0.479        # N/m
        self['pore.thermal_conductivity'] = 8.515   # W/m.K      
        self['pore.vapor_pressure'] = 0.3423        # Pascal
        self['pore.viscosity'] = 1.533E-3           # kg/m.s
                      
if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    air = OpenPNM.Phases.Air(network=pn)