import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM
from OpenPNM.Fluids.__GenericFluid__ import GenericFluid
from OpenPNM.Fluids import models as fm

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
    def __init__(self,name=None,**kwargs):
        super(Water,self).__init__(name=name,**kwargs)
        self._logger.debug("Construct class")
        
    def generate(self):
        self.add_model(propname='pore.diffusivity',
                       model=fm.misc.constant,
                       value=1e-9)
        self.add_model(propname='pore.surface_tension',
                       model=fm.misc.constant,
                       value=0.072)
        self.add_model(propname='pore.contact_angle',
                       model=fm.misc.constant,
                       value=110)
        self.add_model(propname='pore.molar_density',
                       model=fm.misc.constant,
                       value=44445.0)
        self.add_model(propname='pore.viscosity',
                       model=fm.viscosity.reynolds,
                       uo=0.002,
                       b=0.001)
        self.regenerate()  # Include this to allow for old school add_property

if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    water = OpenPNM.Fluids.Water(network=pn)
