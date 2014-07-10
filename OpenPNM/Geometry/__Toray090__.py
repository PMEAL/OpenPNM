"""
module __Toray090__: Subclass of GenericGeometry for a standard Toray TGPH090
gas diffusion layer.s
=============================================================================== 

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)

import OpenPNM
from OpenPNM.Geometry.__GenericGeometry__ import GenericGeometry
from OpenPNM.Geometry import models as gm

class Toray090(GenericGeometry):
    r"""
    Toray090 subclass of GenericGeometry

    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)

    """

    def __init__(self, **kwargs):
        r"""
        Initialize
        """
        super(Toray090,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
        self.add_property(propname='pore.seed',
                          model=gm.pore_seed.random,
                          seed=None)
        self.add_property(propname='throat.seed',
                          model=gm.throat_seed.neighbor,
                          mode='min')
        self.add_property(propname='pore.diameter',
                          model=gm.pore_diameter.sphere,
                          psd_name='weibull_min',
                          psd_shape=2.5,
                          psd_loc=5e-6,
                          psd_scale=4e-6)
        self.add_property(propname='pore.area',
                          model=gm.pore_area.spherical)
        self.add_property(propname='throat.diameter',
                          model=gm.throat_diameter.cylinder,
                          tsd_name='weibull_min',
                          tsd_shape=2.5,
                          tsd_loc=5e-6,
                          tsd_scale=4e-6)                  
        self.add_property(propname='throat.length',
                          model=gm.throat_length.straight)
        self.add_property(propname='throat.volume',
                          model=gm.throat_volume.cylinder)
        self.add_property(propname='throat.area',
                          model=gm.throat_area.cylinder)
        self.add_property(propname='throat.surface_area',
                          model=gm.throat_surface_area.cylinder)
        
if __name__ == '__main__':
    pn = OpenPNM.Network.TestNet()
    test = OpenPNM.Geometry.Stick_and_Ball(loglevel=10,name='test_geom',locations=[0],network=pn)
