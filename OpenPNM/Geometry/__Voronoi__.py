"""
module __Voronoi__: Subclass of GenericGeometry for a standard Geometry created from a Voronoi Diagram
Used with Delaunay Network but could work for others (not tested)
=============================================================================== 

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import OpenPNM.Utilities.transformations as tr
from scipy.spatial import ConvexHull
from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry.__GenericGeometry__ import GenericGeometry
import OpenPNM.Utilities.vertexops as vo

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
    
    def print_throat(self,throats_in):
        r"""
        Print a given throat or list of throats accepted as [1,2,3,...,n]
        e.g geom.print_throat([34,65,99])
        Original vertices plus offset vertices are rotated to align with 
        the z-axis and then printed in 2D
        """
        import matplotlib.pyplot as plt
        throats = []
        for throat in throats_in:
            if throat in range(self.num_throats()):
                throats.append(throat)
            else:
                print("Throat: "+str(throat)+ " not part of geometry")
        if len(throats) > 0:
            verts = self['throat.vertices'][throats]
            offsets = self['throat.offset_vertices'][throats]
            normals = self['throat.normal'][throats]
            coms = self['throat.centroid'][throats]
            for i in range(len(verts)):
                fig = plt.figure()
                vert_2D = tr.rotate_and_chop(verts[i],normals[i],[0,0,1])
                hull = ConvexHull(vert_2D)
                for simplex in hull.simplices:
                    plt.plot(vert_2D[simplex,0], vert_2D[simplex,1], 'k-',linewidth=2)
                plt.scatter(vert_2D[:,0], vert_2D[:,1])
                #centroid = vo.PolyWeightedCentroid2D(vert_2D[hull.vertices])
                offset_2D = tr.rotate_and_chop(offsets[i],normals[i],[0,0,1])
                offset_hull = ConvexHull(offset_2D)
                for simplex in offset_hull.simplices:
                    plt.plot(offset_2D[simplex,0], offset_2D[simplex,1], 'g-',linewidth=2)
                plt.scatter(offset_2D[:,0], offset_2D[:,1])
                #centroid2 = vo.PolyWeightedCentroid2D(offset_2D[offset_hull.vertices])
                " Make sure the plot looks nice by finding the greatest range of points and setting the plot to look square"
                xmax = vert_2D[:,0].max()
                xmin = vert_2D[:,0].min()
                ymax = vert_2D[:,1].max()
                ymin = vert_2D[:,1].min()
                x_range = xmax - xmin
                y_range = ymax - ymin
                if (x_range > y_range):
                    my_range = x_range
                else:
                    my_range = y_range
                lower_bound_x = xmin - my_range*0.5
                upper_bound_x = xmin + my_range*1.5
                lower_bound_y = ymin - my_range*0.5
                upper_bound_y = ymin + my_range*1.5  
                plt.axis((lower_bound_x,upper_bound_x,lower_bound_y,upper_bound_y))
                plt.grid(b=True, which='major', color='b', linestyle='-')
                centroid = tr.rotate_and_chop(coms[i],normals[i],[0,0,1])
                plt.scatter(centroid[0][0],centroid[0][1])
                #plt.scatter(centroid2[0],centroid2[1],c='r')
                fig.show()
        else:
            print("Please provide throat indices")

    def print_pore(self,pores,axis_bounds=[]):
        r"""
        Print all throats around a given pore or list of pores accepted as [1,2,3,...,n]
        e.g geom.print_pore([34,65,99])
        Original vertices plus offset vertices used to create faces and 
        then printed in 3D
        To print all pores (n)
        pore_range = np.arange(0,n-1,1)
        geom.print_pore(pore_range)
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        if len(pores) > 0:
            net_pores = self["pore.map"][pores]
            centroids = self["pore.centroid"][pores]
            #centroids2 = self["pore.com"][pores]
            #for i,pore in enumerate(pores):
            #    centroids[i]=self["pore.centroid"][pore]
            #coords = self._net["pore.coords"][net_pores]
            net_throats = self._net.find_neighbor_throats(pores=net_pores)
            throats = []
            for net_throat in net_throats:
                try:
                    throats.append(self['throat.map'].tolist().index(net_throat))
                except ValueError:
                    " Throat not in this geometry "        
            "Can't create volume from one throat"
            if len(throats)>1:
                verts = self['throat.vertices'][throats]
                normals = self['throat.normal'][throats]
                " Get verts in hull order "
                ordered_verts=[]
                for i in range(len(verts)):
                    vert_2D = tr.rotate_and_chop(verts[i],normals[i],[0,0,1])
                    hull = ConvexHull(vert_2D)
                    ordered_verts.append(verts[i][hull.vertices])
                offsets = self['throat.offset_vertices'][throats]
                "Get domain extents for setting axis "
                if axis_bounds == []:
                    [xmin,xmax,ymin,ymax,zmin,zmax]=vo.vertex_dimension(self._net,pores,parm='minmax')
                else: 
                    [xmin,xmax,ymin,ymax,zmin,zmax]=axis_bounds                
                fig = plt.figure()
                ax = fig.gca(projection='3d')
                outer_items = Poly3DCollection(ordered_verts,linewidths=1, alpha=0.2, zsort='min')
                outer_face_colours=[(1, 0, 0, 0.01)]
                outer_items.set_facecolor(outer_face_colours)
                ax.add_collection(outer_items)
                inner_items = Poly3DCollection(offsets,linewidths=1, alpha=0.2, zsort='min')
                inner_face_colours=[(0, 0, 1, 0.01)]
                inner_items.set_facecolor(inner_face_colours)
                ax.add_collection(inner_items)
                ax.set_xlim(xmin,xmax)
                ax.set_ylim(ymin,ymax)
                ax.set_zlim(zmin,zmax)
                #ax.scatter(coords[:,0],coords[:,1],coords[:,2])
                ax.scatter(centroids[:,0],centroids[:,1],centroids[:,2],c='y')
                #ax.scatter(centroids2[:,0],centroids2[:,1],centroids2[:,2],c='g')
                plt.show()
            else:
                self.print_throat(throats)
        else:
            print("Please provide pore indices")
        
        
if __name__ == '__main__':
    pn = OpenPNM.Network.Delaunay(name='test_net')
    pn.generate(num_pores=100, domain_size=[0.0001,0.0001,0.0001],add_boundaries=True)
    test = OpenPNM.Geometry.Voronoi(loglevel=10,name='test_geom',locations=[0],network=pn)
    test.set_locations(pores=pn.pores('internal'),throats='all') # Don't really know what this does but is needed
    pn.regenerate_geometries()