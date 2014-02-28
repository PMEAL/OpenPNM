"""
module __GenericGeometry__: Base class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)
import OpenPNM

from OpenPNM.Geometry.__GenericGeometry__ import GenericGeometry

class Stick_and_Ball(GenericGeometry):
    r"""
    Stick and Ball subclass of GenericGeometry.

    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)

    """

    def __init__(self, **kwargs):
        r"""
        Initialize
        """
        super(Stick_and_Ball,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
   
        self.add_method(prop='pore_seed',model='random')
        self.add_method(prop='throat_seed',model='neighbor_min')
        self.add_method(prop='pore_diameter',model='sphere',name='weibull_min',shape=2.5,loc=6e-6,scale=2e-5)
        self.add_method(prop='throat_diameter',model='cylinder',name='weibull_min',shape=2.5,loc=6e-6,scale=2e-5)
        self.add_method(prop='pore_volume',model='sphere')
        self.add_method(prop='throat_length',model='straight')
        self.add_method(prop='throat_volume',model='cylinder')
        self.add_method(prop='throat_vector',model='pore_to_pore')
        self.add_method(prop='throat_area',model='cylinder')
        self.add_method(prop='throat_surface_area',model='cylinder')
        
if __name__ == '__main__':
    pn = OpenPNM.Network.TestNet()
    test = OpenPNM.Geometry.Stick_and_Ball(loglevel=10,name='test_geom',locations=[0],network=pn)
