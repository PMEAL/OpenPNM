import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
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
        self.add_property(prop='diffusivity',
                          model='fuller',
                          MA=18,
                          MB=20,
                          vA=1,
                          vB=1)
        self.add_property(prop='molar_density',
                          model='ideal_gas')
        self.add_property(prop='viscosity',
                          model='reynolds',
                          uo=0.0002,
                          b=0.0001)
        self.regenerate()

if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    air = OpenPNM.Fluids.Air(network=pn)