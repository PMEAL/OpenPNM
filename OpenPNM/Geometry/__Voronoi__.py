# -*- coding: utf-8 -*-
"""
===============================================================================
Voronoi --Subclass of GenericGeometry for a standard Geometry created from a
Voronoi Diagram Used with Delaunay Network but could work for others (not tested)
===============================================================================


"""

import scipy as sp
import numpy as np
from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry import GenericGeometry
import time

class Voronoi(GenericGeometry):
    r"""
    Voronoi subclass of GenericGeometry.

    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)

    """

    def __init__(self, fibre_rad = 3e-06,**kwargs):
        r"""
        Initialize
        """
        if int(sp.__version__.split('.')[1]) < 13:
            raise Exception('The installed version of Scipy is too old, Voronoi cannot run')
        super(Voronoi,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
        #if kwargs['fibre_rad']:
        #    fibre_rad = kwargs['fibre_rad']
        #else:
        #    fibre_rad = 3e-06
        self._generate(fibre_rad)

    def _generate(self,fibre_rad):
        r'''
        '''
        timer = []
        elements = []
        timer.append(time.time())
        elements.append("start")
        self.add_model(propname='pore.vertices',
                       model=gm.pore_vertices.voronoi)
        timer.append(time.time())
        elements.append('pore.vertices')
        self.add_model(propname='throat.vertices',
                       model=gm.throat_vertices.voronoi)
        timer.append(time.time())
        elements.append('throat.vertices')
        self.add_model(propname='throat.normal',
                       model=gm.throat_normal.voronoi)
        timer.append(time.time())
        elements.append('throat.normal')
        self.add_model(propname='throat.offset_vertices',
                       model=gm.throat_offset_vertices.distance_transform,
                       offset=fibre_rad,
                       set_dependent = True)
        timer.append(time.time())
        elements.append('throat.offset_vertices')
        
        self._net.trim_occluded_throats()
        
        self.add_model(propname='pore.seed',
                       model=gm.pore_misc.random,
                       seed=self._seed)
        timer.append(time.time())
        elements.append('pore.seed')
        self.add_model(propname='throat.seed',
                       model=gm.throat_misc.neighbor,
                       pore_prop='pore.seed',
                       mode='min')
        timer.append(time.time())
        elements.append('throat.seed')
        self.add_model(propname='pore.volume',
                       model=gm.pore_volume.voronoi)
        timer.append(time.time())
        elements.append('pore.volume')
        self.add_model(propname='pore.diameter',
                       model=gm.pore_diameter.voronoi)
        timer.append(time.time())
        elements.append('pore.diameter')
        self.add_model(propname='pore.area',
                       model=gm.pore_area.spherical)
        timer.append(time.time())
        elements.append('pore.area')
        if 1 == 2:
            self.add_model(propname='throat.area',
                           model=gm.throat_area.voronoi)
            self.add_model(propname='throat.perimeter',
                           model=gm.throat_perimeter.voronoi)
            self.add_model(propname='throat.centroid',
                           model=gm.throat_centroid.centre_of_mass)
            self.add_model(propname='pore.centroid',
                           model=gm.pore_centroid.centre_of_mass)
            self.add_model(propname='pore.indiameter',
                           model=gm.pore_diameter.insphere)
            timer.append(time.time())
            elements.append('pore.indiameter')
            self.add_model(propname='throat.indiameter',
                           model=gm.throat_diameter.incircle) 
            timer.append(time.time())
            elements.append('throat.indiameter')
            
        self.add_model(propname='throat.diameter',
                       model=gm.throat_diameter.voronoi)
        timer.append(time.time())
        elements.append('throat.diameter')
        
        #self.add_model(propname='throat.c2c',
        #               model=gm.throat_length.voronoi)
        #timer.append(time.time())
        #elements.append('throat.c2c')
        self.add_model(propname='throat.length',
                       model=gm.throat_length.constant,
                       const=fibre_rad*2)
        timer.append(time.time())
        elements.append('throat.length')
        self.add_model(propname='throat.volume',
                       model=gm.throat_volume.extrusion)
        timer.append(time.time())
        elements.append('throat.volume')
        self.add_model(propname='throat.surface_area',
                       model=gm.throat_surface_area.extrusion)
        timer.append(time.time())
        elements.append('throat.surface_area')
        timer = sp.asarray(timer)
        print("Geom Timings: ")
        for i in np.arange(1,len(timer)):
            print(elements[i] +" "+ str(timer[i]-timer[i-1]))
        #print(timer[1:]-timer[0:len(timer)-1])
        "Shift the pore coords to the centroids"
        #vo.pore2centroid(self._net)
    
if __name__ == '__main__':
    import OpenPNM
    pn = OpenPNM.Network.Delaunay(name='test_net')
    pn.generate(num_pores=100, domain_size=[0.0001,0.0001,0.0001],add_boundaries=True)
    test = OpenPNM.Geometry.Voronoi(loglevel=10,name='test_geom',locations=[0],network=pn)
    test.set_locations(pores=pn.pores('internal'),throats='all') # Don't really know what this does but is needed
    pn.regenerate_geometries()