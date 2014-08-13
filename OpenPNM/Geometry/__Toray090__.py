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
from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry.__GenericGeometry__ import GenericGeometry

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
        self._generate()
                       
    def _generate(self):
        r'''
        '''        
        self.add_model(propname='pore.seed',
                       model=gm.pore_misc.random,
                       seed=None)
        self.add_model(propname='throat.seed',
                       model=gm.throat_misc.neighbor,
                       pore_prop='pore.seed',
                       mode='min')
        self.add_model(propname='pore.diameter',
                       model=gm.pore_diameter.sphere_from_radius,
                       psd_name='weibull_min',
                       psd_shape=1.5,
                       psd_loc=14e-6,
                       psd_scale=2e-6)
        self.add_model(propname='pore.area',
                       model=gm.pore_area.spherical)
        self.add_model(propname='pore.volume',
                       model=gm.pore_volume.sphere)
        self.add_model(propname='throat.diameter',
                       model=gm.throat_diameter.cylinder_from_radius,
                       tsd_name='weibull_min',
                       tsd_shape=1.5,
                       tsd_loc=14e-6,
                       tsd_scale=2e-6)                  
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
