"""
module __SGL10__: Subclass of GenericGeometry for an SGL10 gas diffusion 
layer
===============================================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)
import OpenPNM

from OpenPNM.Geometry.__GenericGeometry__ import GenericGeometry

class SGL10(GenericGeometry):
    r"""
    SGL10 subclass of GenericGeometry.

    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=Info, 30=Warning, 40=Error, 50=Critical)

    """

    def __init__(self, **kwargs):
        r"""
        Initialize
        """
        super(SGL10,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
   
        self.add_property(prop='pore_seed',model='random')
        self.add_property(prop='throat_seed',model='neighbor',mode='min')
        self.add_property(prop='pore_diameter',model='sphere',psd_name='weibull_min',psd_shape=2.5,psd_loc=9e-6,psd_scale=5e-6)
        self.add_property(prop='throat_diameter',model='cylinder',tsd_name='weibull_min',tsd_shape=2.5,tsd_loc=9e-6,tsd_scale=5e-5)
        self.add_property(prop='pore_volume',model='sphere')
        self.add_property(prop='throat_length',model='straight')
        self.add_property(prop='throat_volume',model='cylinder')
        self.add_property(prop='throat_vector',model='pore_to_pore')
        self.add_property(prop='throat_area',model='cylinder')
        self.add_property(prop='throat_surface_area',model='cylinder')
        
if __name__ == '__main__':
    pn = OpenPNM.Network.TestNet()
    test = OpenPNM.Geometry.Stick_and_Ball(loglevel=10,name='test_geom',locations=[0],network=pn)
