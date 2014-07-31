"""
module __StickBall__: Subclass of GenericGeometry for a standard 'stick & ball'
geometrical model
==================================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)

import OpenPNM
from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry.__GenericGeometry__ import GenericGeometry

class Stick_and_Ball(GenericGeometry):
    r"""
    Stick and Ball subclass of GenericGeometry.  This subclass is meant as a 
    basic default geometry to get started quickly.  

    Parameters
    ----------
    name : string
        The name of the object, which is also used as the label where this 
        geometry is defined.

    """

    def __init__(self, **kwargs):
        r"""
        Initialize
        """
        super(Stick_and_Ball,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
        
    def generate(self):   
        self.add_model(propname='pore.seed',
                       model=gm.pore_seed.random,
                       seed=None)
        self.add_model(propname='throat.seed',
                       model=gm.throat_seed.neighbor,
                       mode='min')
        self.add_model(propname='pore.diameter',
                       model=gm.pore_diameter.sphere,
                       psd_name='weibull_min',
                       psd_shape=2.5,
                       psd_loc=0,
                       psd_scale=0.5)
        self.add_model(propname='pore.area',
                       model=gm.pore_area.spherical)
        self.add_model(propname='pore.volume',
                       model=gm.pore_volume.sphere)
        self.add_model(propname='throat.diameter',
                       model=gm.throat_diameter.cylinder,
                       tsd_name='weibull_min',
                       tsd_shape=2.5,
                       tsd_loc=0,
                       tsd_scale=0.5)                  
        self.add_model(propname='throat.length',
                       model=gm.throat_length.straight)
        self.add_model(propname='throat.volume',
                       model=gm.throat_volume.cylinder)
        self.add_model(propname='throat.area',
                       model=gm.throat_area.cylinder)
        self.add_model(propname='throat.surface_area',
                       model=gm.throat_surface_area.cylinder)
        
if __name__ == '__main__':
    pn = OpenPNM.Network.TestNet()
    pass
