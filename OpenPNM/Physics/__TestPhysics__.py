"""
module Physics
===============================================================================

"""

import OpenPNM
import scipy as sp
from OpenPNM.Physics import models as pm
from OpenPNM.Physics.__GenericPhysics__ import GenericPhysics

class TestPhysics(GenericPhysics):
    r"""
    Base class to generate a generic Physics object.  The user must specify models
    and parameters for the all the properties they require. Classes for several
    common Physics are included with OpenPNM and can be found under OpenPNM.Physics.

    Parameters
    ----------
    network : OpenPNM Network object 
        The network to which this Physics should be attached
        
    phase : OpenPNM Phase object 
        The Phase object to which this Physics applies
        
    pores and throats : array_like
        The pores and throats where this Physics object applies

    """

    def __init__(self,**kwargs):
        super(TestPhysics,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        self._generate()
        
    def _generate(self):
        for phase in self._phases:
            temp = [item.split('.')[1] for item in phase.props()]
            if 'viscosity' in temp:
                self['throat.hydraulic_conductance'] = 1
            if 'diffusivity' in temp:
                self['throat.diffusive_conductance'] = 1
            if 'surface_tension' in temp:
                self['throat.capillary_pressure'] = 1/self._Phases[0]['throat.diameter']
        
if __name__ == '__main__':
    print('none yet')


