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
        self.add_property(prop='throat_seed',model='neighbor_min')
        self.add_property(prop='pore_diameter',model='sphere',name='weibull_min',shape=.987,loc=1.58e-5,scale=5.01e-6)
        self.add_property(prop='throat_diameter',model='cylinder',name='weibull_min',shape=.987,loc=1.58e-5,scale=5.01e-6)
        self.add_property(prop='pore_volume',model='cube')
        self.add_property(prop='throat_length',model='straight')
        self.add_property(prop='throat_volume',model='cuboid')
        self.add_property(prop='throat_vector',model='pore_to_pore')
        self.add_property(prop='throat_area',model='cuboid')
        self.add_property(prop='throat_surface_area',model='cylinder')
        
if __name__ == '__main__':
    pn = OpenPNM.Network.TestNet()
    test = OpenPNM.Geometry.Stick_and_Ball(loglevel=10,name='test_geom',locations=[0],network=pn)
