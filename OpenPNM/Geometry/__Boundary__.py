"""
module __Boundary__: Subclass of GenericGeometry for Boundary Pores
==================================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)
import OpenPNM

from OpenPNM.Geometry.__GenericGeometry__ import GenericGeometry

class Boundary(GenericGeometry):
    r"""
    Boundary subclass of GenericGeometry.

    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)

    """

    def __init__(self, **kwargs):
        r"""
        Initialize
        """
        super(Boundary,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
   
        self.add_property(prop='pore_seed',model='constant',value=1.0)
        self.add_property(prop='throat_seed',model='constant',value=1.0)
        self.add_property(prop='pore_diameter',model='constant',value=0)
        self.add_property(prop='throat_diameter',model='constant',value=0)
        self.add_property(prop='pore_volume',model='constant',value=0.0)
        self.add_property(prop='throat_length',model='constant',value=0.0)
        self.add_property(prop='throat_volume',model='constant',value=0.0)
        self.add_property(prop='throat_vector',model='pore_to_pore')
        self.add_property(prop='throat_area',model='cylinder')
        self.add_property(prop='throat_surface_area',model='constant',value=0.0)
        
if __name__ == '__main__':
    pn = OpenPNM.Network.TestNet()
    test = OpenPNM.Geometry.Boundary(loglevel=10,name='test_geom',locations=[0],network=pn)
