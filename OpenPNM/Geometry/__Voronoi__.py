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
import scipy as sp
import numpy as np
import _transformations as tr
from scipy.spatial import ConvexHull
from math import atan2

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
        self._add_throat_props() # This sets the key throat data for calculating pore and throat properties later
        self.add_property(prop='pore_seed',model='random')
        self.add_property(prop='throat_seed',model='neighbor_min')
        self.add_property(prop='pore_volume',model='voronoi2') # Volume must come before diameter
        self.add_property(prop='pore_diameter',model='voronoi')
        self.add_property(prop='pore_centroid',model='voronoi')
        # Throat must come before Pore to get the offset vertices
        #self.add_property(prop='throat_area', model='voronoi') # Area must come before diameter, also sets face centroid and perimeter
        self.add_property(prop='throat_diameter',model='voronoi')
        self.add_property(prop='throat_centroid',model='voronoi')
        #self.add_property(prop='throat_length',model='voronoi') # Length must come before volume
        self.add_property(prop='throat_length',model='constant',value=1e-06)
        self.add_property(prop='throat_volume',model='voronoi')
        self.add_property(prop='throat_vector',model='pore_to_pore') # Not sure how to do this for centre to centre as we might need to split into two vectors
        self.add_property(prop='throat_surface_area',model='voronoi')

    def _add_throat_props(self,radius=1e-06):
        r"""
        This method does all the throat properties for the voronoi cages 
        including offseting the vertices surrounding each pore by an amount that
        replicates erroding the facet of each throat by the fibre radius 

        For each connection or throat find the shared vertices
        Rotate the vertices to align with the xy-plane and get rid of z-coordinate
        Merge vertices within a certain distance of each other based on the range of data points
        Compute the convex hull of the 2D points giving a set of simplices which define neighbouring vertices in a clockwise fashion
        For each triplet calculate the offset position given the fibre radius
        Translate back into 3D
        """
        connections = self._net['throat.conns']
        coords = self._net['pore.coords']
        verts = self._net['pore.vertices']
        normals = coords[connections[:,0]]-coords[connections[:,1]]
        area = sp.ndarray(len(connections),dtype=object)
        perimeter = sp.ndarray(len(connections),dtype=object)
        #centroids = sp.ndarray(len(connections),dtype=object)
        offset_verts = sp.ndarray(len(connections),dtype=object)
        for i,throat_pair in enumerate(connections):
            shared_verts = []
            pore_a = throat_pair[0]
            pore_b = throat_pair[1]
            " Identify shared verts "
            for vert_a in verts[pore_a]:
                for vert_b in verts[pore_b]:
                    if (vert_a[0] == vert_b[0]) and (vert_a[1] == vert_b[1]) and (vert_a[2] == vert_b[2]):
                        shared_verts.append(vert_a)
            if len(shared_verts) >=3:
                shared_verts = np.asarray(shared_verts)
                area[i],perimeter[i],offset_verts[i] = self._get_throat_geom(shared_verts,normals[i])
                #centroids[i]=self._centroid(shared_verts)
            else:
                area[i]=0.0

        self._net.set_data(prop='area',throats='all',data=area)
        #self._net['throat.area']=area
        #self._net.set_data(prop='perimeter',throats='all',data=perimeter)
        self._net['throat.perimeter']=perimeter
        #self.set_throat_data(prop='centroid',data=centroids)
        #self._net.set_throat_data(prop='offset_verts',data=offset_verts)
        self._net['throat.offset_verts']=offset_verts
        " Temporary Code to remove the smallest 2% of throat area connections "
        " This can be used in algorithms to ignore certain connections if required - like in the range of capillary pressures in OP "
        #smallest_throat = min(area)
        #largest_throat = max(area)
        average_area = sp.mean(area)
        cutoff = average_area/100
        #remove_range = smallest_throat + ((largest_throat-smallest_throat)*0.02)
        excluded_throats = []
        for i,throat in enumerate(connections):
            if area[i]<=cutoff:
            #if area[i]==0.0:
                excluded_throats.append(i)
        excluded_throats = np.asarray(excluded_throats)
        if len(excluded_throats) > 0:
            self._net.trim(throats=excluded_throats)
    
    def _get_throat_geom(self,verts,normal):
        z_axis = [0,0,1]
        " For boundaries some facets will already be aligned with the axis - if this is the case a rotation is unnecessary and could also cause problems "
        angle = tr.angle_between_vectors(normal,z_axis)
        if (angle==0.0)or(angle==np.pi):
            "We are already aligned"
            rotate_input = False
            facet = verts
        else:
            rotate_input = True
            M = tr.rotation_matrix(tr.angle_between_vectors(normal,z_axis),tr.vector_product(normal,z_axis))
            facet = np.dot(verts,M[:3,:3].T)
        x = facet[:,0]
        y = facet[:,1]
        z = facet[:,2]
        if (np.around(z.std(),3)!=0.000):
            print("Rotation failed")
        facet_coords_2D = np.column_stack((x,y))
        hull = ConvexHull(facet_coords_2D)
        verts_2D = facet_coords_2D[hull.vertices]
        fused_verts=self._fuse_verts(verts_2D)
        if len(fused_verts) <3:
            #we fused too many
            fused_verts=self._fuse_verts(verts_2D,0.0025)
        offset = []
        fibre_rad = 3e-06
        for i,vert in enumerate(fused_verts):
            " Collect three adjacent points and compute the offset of the first "
            triplet = (vert, np.roll(fused_verts,-1,axis=0)[i],np.roll(fused_verts,1,axis=0)[i])
            offset.append(self._offset_vertex(triplet,fibre_rad))
        offset = np.asarray(offset)
    
        original_area = self._PolyArea2D(verts_2D)            
        all_points = np.concatenate((verts_2D,offset),axis=0)
        try:
            total_hull = ConvexHull(all_points,qhull_options='Pp') #ignores very small angles
            total_area = self._PolyArea2D(all_points[total_hull.vertices])
        except sp.spatial.qhull.QhullError:
            print(all_points)
            total_area =999

        if (total_area>original_area): # Throat is fully occluded
            throat_area = 0.0
            throat_perimeter = 0.0
            output_offset = []
        else:
            offset_hull = ConvexHull(offset)
            offset_verts_2D = offset[offset_hull.vertices]
            throat_area = self._PolyArea2D(offset_verts_2D)
            throat_perimeter = self._PolyPerimeter2D(offset_verts_2D)
            " Make 3D again in rotated plane "
            offset_verts_3D = np.column_stack((offset_verts_2D,z[0:len(offset_verts_2D)]))
            " Get matrix to un-rotate the co-ordinates back to the original orientation if we rotated in the first place"
            if (rotate_input):
                M1 = tr.inverse_matrix(M)
                " Unrotate the offset coordinates "
                output_offset = np.dot(offset_verts_3D,M1[:3,:3].T)
            else:
                output_offset = offset_verts_3D

        return throat_area, throat_perimeter, output_offset
    
    def _offset_vertex(self,points,rad = 0.01):
        " We are passed in a set of 3 points forming vertices of two adjoining simplexes of the convex hull of a voronoi facet "
        " We need to offset the vertices normal to the fibre direction (or adjoining vectors) by the fibre radius "
        " This is achieved by finding the half angle between the two adjoining vectors and a direction "
        " Mid-point must be the first in the array "
        p0 = np.array(points[0])
        p1 = np.array(points[1])
        p2 = np.array(points[2])
        " Now make the midpoint the origin "
        vector1 = p1-p0
        vector2 = p2-p0

        " Find what quadrant the vector is pointing in - atan2 function takes account of signs "
        " 0 means aligned with x-axis, pi is aligned with -xaxis, positive numbers are positive y and negative numbers are negative y "
        " The angle between the vectors should always be within 180 degrees of one another in a convex hull "

        q1 = atan2(vector1[1],vector1[0])
        q2 = atan2(vector2[1],vector2[0])
        alpha = 0.5*tr.angle_between_vectors(vector1,vector2)
    
        " We always want to offset from the first vertex we get to - going anti-clockwise from the x-axis "
        " Check if both vectors point up or both point down - if so the first one we get to will have smaller q value "
        if q1*q2 >=0.0:
            if q1<q2:
                theta = q1
            else:
                theta = q2
        else:
            "if vector 1 is more rotated away from positive xaxis than vector 2 is rotated away from negative xaxis - use it"
            " and vice-versa "
            if (abs(q1)+abs(q2)>np.pi):
                "vectors are pointing negative x so take whichever has positive q-value - like a pacman facing left"
                if q1>=0:
                    theta = q1
                else:
                    theta = q2
            else:
                "vectors are pointing positive x so take whichever is negative"
                if q1<=0:
                    theta = q1
                else:
                    theta = q2
    
        x = rad*np.cos(alpha+theta)/np.sin(alpha)
        y = rad*np.sin(alpha+theta)/np.sin(alpha)

        "Add the midpoint back in"
        output = [x+p0[0],y+p0[1]]    
    
        return output

    def _dist2(self,p1, p2):
        return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2

    def _fuse(self,points, d):
        ret = []
        d2 = d * d
        n = len(points)
        taken = [False] * n
        for i in range(n):
            if not taken[i]:
                count = 1
                point = [points[i][0], points[i][1]]
                taken[i] = True
                for j in range(i+1, n):
                    if self._dist2(points[i], points[j]) < d2:
                        point[0] += points[j][0]
                        point[1] += points[j][1]
                        count+=1
                        taken[j] = True
                point[0] /= count
                point[1] /= count
                ret.append((point[0], point[1]))
        return ret

    def _fuse_verts(self,verts,percentage=0.05):
        #Work out largest span
        x_span = max(verts[:,0])- min(verts[:,0])
        y_span = max(verts[:,1])- min(verts[:,1])
        if x_span > y_span:
            tolerance = x_span*percentage
        else:
            tolerance = y_span*percentage
        #fuse vertices lying within 5% of the largest span    
        return self._fuse(verts,tolerance)

    def _PolyArea2D(self,pts):
        lines = np.hstack([pts,np.roll(pts,-1,axis=0)])
        area = 0.5*abs(sum(x1*y2-x2*y1 for x1,y1,x2,y2 in lines))
        return area

    def _PolyPerimeter2D(self,pts):
        lines = np.hstack([pts,np.roll(pts,-1,axis=0)])
        perimeter = sum(np.sqrt((x2-x1)**2+(y2-y1)**2) for x1,y1,x2,y2 in lines)
        return perimeter
        
if __name__ == '__main__':
    pn = OpenPNM.Network.Delaunay(name='test_net')
    pn.generate(num_pores=100, domain_size=[0.0001,0.0001,0.0001],add_boundaries=True)
    test = OpenPNM.Geometry.Voronoi(loglevel=10,name='test_geom',locations=[0],network=pn)
    test.set_locations(pores=pn.pores('internal'),throats='all') # Don't really know what this does but is needed
    pn.regenerate_geometries()