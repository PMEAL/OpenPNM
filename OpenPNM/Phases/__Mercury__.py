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
        self['pore.density'] = 7000  # kg/m3
        self['pore.contact_angle'] = 140  # Degrees
        self['pore.surface_tension'] = 0.480  # N/m

if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    air = OpenPNM.Phases.Air(network=pn)