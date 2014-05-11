import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM
from OpenPNM.Fluids.__GenericFluid__ import GenericFluid

class Water(GenericFluid):
    r'''
    Creates Fluid object with a default name 'water' and preset values for water
    
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
    def __init__(self,name='water',**kwargs):
        super(Water,self).__init__(name=name,**kwargs)
        self._logger.debug("Construct class")
        self.set_pore_data(prop='Tc',data=647.096)
        self.set_pore_data(prop='Pc',data=22.06e6)
        self.set_pore_data(prop='MW',data=0.0291)
        self.add_property(prop='diffusivity',model='constant',value=2e-9)
        self.add_property(prop='viscosity',model='constant',value=0.001)
        self.add_property(prop='molar_density',model='constant',value=44445)
        self.add_property(prop='surface_tension',model='constant',value=0.072)
        self.add_property(prop='contact_angle',model='constant',value=110)
        self.regenerate()

if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    water = OpenPNM.Fluids.Water(network=pn)
