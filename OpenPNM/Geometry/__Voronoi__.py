"""
module __Voronoi__: Subclass of GenericGeometry for a standard Geometry created from a Voronoi Diagram
Used with Delaunay Network but could work for others (not tested)
=============================================================================== 

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import OpenPNM
import scipy as sp
from OpenPNM.Geometry import models as gm
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
        if int(sp.__version__.split('.')[1]) < 13:
            raise Exception('The installed version of Scipy is too old, Voronoi cannot run')
        super(Voronoi,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
        if kwargs['fibre_rad']:
            fibre_rad = kwargs['fibre_rad']
        else:
            fibre_rad = 3e-06
        self._generate(fibre_rad)
    
    def _generate(self,fibre_rad):
        r'''
        ''' 
        self.add_model(propname='pore.vertices',
                       model=gm.pore_vertices.voronoi)
        self.add_model(propname='throat.vertices',
                       model=gm.throat_vertices.voronoi)
        self.add_model(propname='throat.normal',
                       model=gm.throat_normal.voronoi)
        self.add_model(propname='throat.offset_vertices',
                       model=gm.throat_offset_vertices.voronoi,
                       offset=fibre_rad)
        self.add_model(propname='pore.seed',
                       model=gm.pore_misc.random,
                       seed=self._seed)
        self.add_model(propname='throat.seed',
                       model=gm.throat_misc.neighbor,
                       pore_prop='pore.seed',
                       mode='min')
        self.add_model(propname='pore.volume',
                       model=gm.pore_volume.voronoi)
        self.add_model(propname='pore.diameter',
                       model=gm.pore_diameter.voronoi)
        #self.add_model(propname='pore.centroid',
        #               model=gm.pore_centroid.voronoi)
        self.add_model(propname='pore.area',
                       model=gm.pore_area.spherical)
        self.add_model(propname='throat.area',
                       model=gm.throat_area.voronoi)
        #self._net.trim_occluded_throats()
        self.add_model(propname='throat.perimeter',
                       model=gm.throat_perimeter.voronoi)
        self.add_model(propname='throat.shape_factor',
                       model=gm.throat_shape_factor.compactness)
        #self.add_model(propname='throat.centroid',
        #               model=gm.throat_centroid.voronoi)
        self.add_model(propname='throat.centroid',
                       model=gm.throat_centroid.centre_of_mass)
        self.add_model(propname='pore.centroid',
                       model=gm.pore_centroid.centre_of_mass)
        self.add_model(propname='pore.diameter',
                       model=gm.pore_diameter.voronoi)
        self.add_model(propname='pore.indiameter',
                       model=gm.pore_diameter.insphere)
        self.add_model(propname='throat.diameter',
                       model=gm.throat_diameter.voronoi)
        self.add_model(propname='throat.indiameter',
                       model=gm.throat_diameter.incircle) 
        self.add_model(propname='throat.c2c',
                       model=gm.throat_length.voronoi)
        self.add_model(propname='throat.length',
                       model=gm.throat_length.constant,
                       const=fibre_rad*2)
        self.add_model(propname='throat.volume',
                       model=gm.throat_volume.extrusion)
        self.add_model(propname='throat.surface_area',
                       model=gm.throat_surface_area.extrusion)
        "Shift the pore coords to the centroids"
        #vo.pore2centroid(self._net)
    
if __name__ == '__main__':
    pn = OpenPNM.Network.Delaunay(name='test_net')
    pn.generate(num_pores=100, domain_size=[0.0001,0.0001,0.0001],add_boundaries=True)
    test = OpenPNM.Geometry.Voronoi(loglevel=10,name='test_geom',locations=[0],network=pn)
    test.set_locations(pores=pn.pores('internal'),throats='all') # Don't really know what this does but is needed
    pn.regenerate_geometries()