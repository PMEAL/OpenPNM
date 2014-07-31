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
        self._generate()
        
    def _generate(self):
        self['pore.diffusivity'] = 1e-9
        self['pore.surface_tension'] = 0.072
        self['pore.contact_angle'] = 110.0
        self['pore.molar_density'] = 44445.0
        self.add_model(propname='pore.viscosity',
                       model=fm.viscosity.reynolds,
                       uo=0.002,
                       b=0.001)

if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    water = OpenPNM.Fluids.Water(network=pn)
