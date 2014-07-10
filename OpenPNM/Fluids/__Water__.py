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
        self.add_property(propname='pore.diffusivity',
                          model=fm.diffusivity.fuller,
                          MA=18,
                          MB=20,
                          vA=1,
                          vB=1)
        self.add_property(propname='pore.surface_tension',
                          model=fm.misc.constant,
                          value=0.072)
        self.add_property(propname='pore.contact_angle',
                          model=fm.misc.constant,
                          value=110)
        self.add_property(propname='pore.molar_density',
                          model=fm.molar_density.ideal_gas,
                          sigma_sl=0.02,
                          sigma_sg=0.09)
        self.add_property(propname='pore.viscosity',
                          model=fm.viscosity.reynolds,
                          uo=0.002,
                          b=0.001)

if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    water = OpenPNM.Fluids.Water(network=pn)
