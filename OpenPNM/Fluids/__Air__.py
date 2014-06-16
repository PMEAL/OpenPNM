import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM
from OpenPNM.Fluids.__GenericFluid__ import GenericFluid

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
        self.set_pore_data(prop='Tc',data=132.65)
        self.set_pore_data(prop='Pc',data=3.771e6)
        self.set_pore_data(prop='MW',data=0.0291)
        self.add_property(prop='diffusivity',model='Fuller',MA=0.03199,MB=0.0291,vA=16.3,vB=19.7)
        self.add_property(prop='viscosity',model='constant',value=1.9e-5)
        self.add_property(prop='molar_density',model='ideal_gas',R=8.314)
        self.add_property(prop='surface_tension',model='constant',value=0.072)
        self.add_property(prop='contact_angle',model='constant',value=70)
        self.regenerate()

if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    air = OpenPNM.Fluids.Air(network=pn)