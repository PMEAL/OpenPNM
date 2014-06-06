"""
module __Voronoi__: Subclass of GenericGeometry for a standard Geometry created from a Voronoi Diagram
Used with Delaunay Network but could work for others (not tested)
=============================================================================== 

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)
import OpenPNM

from OpenPNM.Geometry.__GenericGeometry__ import GenericGeometry

class Voronoi(GenericGeometry):
    r"""
    Voronoi subclass of GenericGeometry.

    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)

    """

    def __init__(self, **kwargs):
        r"""
        Initialize
        """
        super(Voronoi,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
   
        self.add_property(prop='pore_seed',model='random')
        self.add_property(prop='throat_seed',model='neighbor_min')
        self.add_property(prop='pore_volume',model='voronoi2') # Volume must come before diameter
        self.add_property(prop='pore_diameter',model='voronoi')
        # Throat must come before Pore to get the offset vertices
        #self.add_property(prop='throat_area', model='voronoi') # Area must come before diameter, also sets face centroid and perimeter
        self.add_property(prop='throat_diameter',model='voronoi')
        #self.add_property(prop='throat_length',model='voronoi') # Length must come before volume
        self.add_property(prop='throat_length',model='constant',value=1e-06)
        self.add_property(prop='throat_volume',model='voronoi')
        self.add_property(prop='throat_vector',model='pore_to_pore') # Not sure how to do this for centre to centre as we might need to split into two vectors
        self.add_property(prop='throat_surface_area',model='voronoi')

        
if __name__ == '__main__':
    pn = OpenPNM.Network.Delaunay(name='test_net')
    pn.generate(num_pores=100, domain_size=[0.0001,0.0001,0.0001],add_boundaries=True)
    test = OpenPNM.Geometry.Voronoi(loglevel=10,name='test_geom',locations=[0],network=pn)
    test.set_locations(pores=pn.pores('internal'),throats='all') # Don't really know what this does but is needed
    pn.regenerate_geometries()