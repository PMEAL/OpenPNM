"""
module Physics
===============================================================================

"""
import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)

import OpenPNM
import scipy as sp
from OpenPNM.Physics import models as pm
from OpenPNM.Physics.__GenericPhysics__ import GenericPhysics

class Standard(GenericPhysics):
    r"""
    Base class to generate a generic Physics object.  The user must specify models
    and parameters for the all the properties they require. Classes for several
    common Physics are included with OpenPNM and can be found under OpenPNM.Physics.

    Parameters
    ----------
    network : OpenPNM Network object 
        The network to which this Physics should be attached
        
    fluid : OpenPNM Fluid object 
        The Fluid object to which this Physics applies
        
    geometry : OpenPNM Geometry object
        The Geometry object to which this Physics applies
    
    name : str
        A unique string name to identify the Physics object, typically same as 
        instance name but can be anything.
    
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    
    loggername : string (optional)
        Sets a custom name for the logger, to help identify logger messages

    """

    def __init__(self,network,fluid,geometry,name=None,**kwargs):
        super(Standard,self).__init__(network=network,fluid=fluid,geometry=geometry,name=name)
        self._logger.debug("Construct class")
        temp = [item.split('.')[1] for item in fluid.props()]
        if 'viscosity' in temp:
            self.add_property(propname='throat.hydraulic_conductance',
                              model=pm.hydraulic_conductance.hagen_poiseuille)
        if 'diffusivity' in temp:
            self.add_property(propname='throat.diffusive_conductance',
                              model=pm.diffusive_conductance.bulk_diffusion)
        if 'surface_tension' in temp:
            self.add_property(propname='throat.capillary_pressure',
                              model=pm.capillary_pressure.washburn)

        







if __name__ == '__main__':
    print('none yet')


