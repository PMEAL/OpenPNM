import sys, os
import OpenPNM
from OpenPNM.Fluids.__GenericFluid__ import GenericFluid
from OpenPNM.Fluids import models as fm

class Air(GenericFluid):
    r"""
    Creates Fluid object with a default name 'air' and preset values for air
    
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
    >>> air = OpenPNM.Fluids.Air(network=pn)
    """
    def __init__(self,name=None,**kwargs):
        super(Air,self).__init__(name=name,**kwargs)
        self._logger.debug("Construct class")
        self._generate()
        
    def _generate(self):
        self.add_model(propname='pore.diffusivity',
                       model=fm.diffusivity.fuller,
                       MA=18,
                       MB=20,
                       vA=1,
                       vB=1)
        self.add_model(propname='pore.molar_density',
                       model=fm.molar_density.ideal_gas)
        self.add_model(propname='pore.viscosity',
                       model=fm.viscosity.reynolds,
                       uo=9.16656E-6,
                       b=-2.34621E-3)
        self['pore.molecular_weight'] = 28.97
        self['pore.critical_temperature'] = 132.5
        self['pore.critical_volume'] = 0.0883

if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    air = OpenPNM.Fluids.Air(network=pn)