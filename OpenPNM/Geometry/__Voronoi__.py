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
from scipy.spatial import Delaunay

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
        self._add_throat_props(radius=2e-06) # This sets the key throat data for calculating pore and throat properties later
        self.add_property(prop='pore_seed',model='random')
        self.add_property(prop='throat_seed',model='neighbor_min')
        self.add_property(prop='pore_volume',model='voronoi') # Volume must come before diameter
        self.add_property(prop='pore_diameter',model='voronoi')
        self.add_property(prop='pore_centroid',model='voronoi')
        self.add_property(prop='throat_diameter',model='voronoi')
        self.add_property(prop='throat_centroid',model='voronoi')
        self.add_property(prop='throat_length',model='constant',value=2e-06)
        self.add_property(prop='throat_volume',model='voronoi')
        self.add_property(prop='throat_vector',model='pore_to_pore') # Not sure how to do this for centre to centre as we might need to split into two vectors
        self.add_property(prop='throat_surface_area',model='voronoi')

    def _add_throat_props(self,radius=1e-06):
        r"""
        Main Loop         
        This method does all the throat properties for the voronoi cages 
        including calling the offseting routine which offsets the vertices surrounding each pore by an amount that
        replicates erroding the facet of each throat by the fibre radius 
        """
        connections = self._net['throat.conns']
        coords = self._net['pore.coords']
        verts = self._net['pore.vertices']
        normals = coords[connections[:,0]]-coords[connections[:,1]]
        area = sp.ndarray(len(connections),dtype=object)
        perimeter = sp.ndarray(len(connections),dtype=object)
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
                area[i],perimeter[i],offset_verts[i] = self._get_throat_geom(shared_verts,normals[i],radius)
            else:
                area[i]=0.0

        self._net.set_data(prop='area',throats='all',data=area)
        self._net['throat.perimeter']=perimeter
        self._net['throat.offset_verts']=offset_verts
        " Temporary Code to remove the smallest 2% of throat area connections "
        " This can be used in algorithms to ignore certain connections if required - like in the range of capillary pressures in OP "

        average_area = sp.mean(area)
        cutoff = average_area/100
        excluded_throats = []
        for i,throat in enumerate(connections):
            if area[i]<=cutoff:
                excluded_throats.append(i)
        excluded_throats = np.asarray(excluded_throats)
        if len(excluded_throats) > 0:
            self._net.trim(throats=excluded_throats)
    
    def _get_throat_geom(self,verts,normal,fibre_rad):
        r"""
        For one set of vertices defining a throat return the key properties
        This is the main loop for calling other sub-routines.
        General Method:
            For each connection or throat defined by the shared vertices
            Rotate the vertices to align with the xy-plane and get rid of z-coordinate
            Compute the convex hull of the 2D points giving a set of simplices which define neighbouring vertices in a clockwise fashion
            For each triplet calculate the offset position given the fibre radius
            Check for overlapping vertices and ones that lie outside the original hull - recalculate position to offset from or ignore if all overlapping
            Calculate Area and Perimeter if successfully generated offset vertices to replicate eroded throat            
            Translate back into 3D
        Any Errors encountered result in the throat area being zero and no vertices being passed back
        These Errors are not coding mistakes but failures to obtain an eroded facet with non-zero area:
        Error 1: Less than 3 vertices in the convex hull - Should never happen (unless 2 points are incredibly close together)
        Error 2: The largest span of points is less than twice the fibre radius (i.e. throat will definitley be occluded)
        Error 3: All the offset vertices overlap with at least one other vertex - Throat fully occluded
        Error 4: Not enough offset vertices to continue - Throat fully occluded
        Error 5: An offset vertex is outside the original set of points - Throat fully occluded
        """        
        z_axis = [0,0,1]
        throat_area = 0.0
        throat_perimeter = 0.0
        output_offset = []
        Error = 0
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
        " Work out span of points and set axes scales to cover this and be equal in both dimensions "
        x_range = x.max() - x.min()
        y_range = y.max() - y.min()
        if (x_range > y_range):
            my_range = x_range
        else:
            my_range = y_range
        if (np.around(z.std(),3)!=0.000):
            print("Rotation failed")
        facet_coords_2D = np.column_stack((x,y))
        hull = ConvexHull(facet_coords_2D)
        verts_2D = facet_coords_2D[hull.vertices]
        offset = self._outer_offset(verts_2D,fibre_rad)
        " At this point we may have overlapping areas for which we need to offset from a new point "
        overlap_array,sweep_radius,line_points = self._set_overlap(verts_2D,offset)
        temp_vert_list=[]
        #new_vert_list=[]
        if (len(verts_2D) <3):
            "Error: Fused Too Many Verts"
            Error = 1
        elif(my_range < fibre_rad*2):
            "Error: Facet Too small to Erode"
            Error = 2
        else:
            if overlap_array.any()==False:
                " If no overlaps don't worry"
                "no overlaps"
            elif self._all_overlap(overlap_array)==True:
                " If all overlaps then throat is fully occluded"
                "Error: Throat fully occluded"
                Error = 3
            else:
                " If one or two sets of overlaps exist and at least one vertex is not overlapped then we need to do a bit more work "
                " Do some linalg to find a new point to offset from saving un-overlapped verts and newly created verts in a temporary list "
 
                for i in range(np.shape(line_points)[0]):
                    if np.sum(overlap_array[i])==0.0:
                        temp_vert_list.append(verts_2D[i])
                    else:
                        my_lines=[]
                        for j in range(np.shape(line_points)[0]):
                        
                            if overlap_array[i][j] ==1 and overlap_array[j][i]==1:
                                list_a = line_points[i][j]
                                list_b = line_points[j][i]
                                my_lines = self._symmetric_difference(list_a,list_b)
 
                        my_lines=np.asarray(my_lines)
 
                        if len(my_lines)==2:
                            try:
                                quad_points=verts_2D[my_lines]
                            except IndexError:
                                print(verts_2D, my_lines)
                            my_new_point = self._new_point(quad_points)
                            temp_vert_list.append(my_new_point)
                            #new_vert_list.append(my_new_point)
                        
                verts_2D=np.asarray(self._unique_list(temp_vert_list))
                #new_vert_list=np.asarray(self._unique_list(new_vert_list))
                if len(verts_2D) >=3:
                    offset = self._outer_offset(verts_2D,fibre_rad)
                    overlap_array,sweep_radius,line_points = self._set_overlap(verts_2D,offset)
                else:
                    Error = 4

        if len(offset) >= 3 and Error == 0:    
            " Now also check whether any of the offset points lie outside the original convex hull "
            original_area = self._PolyArea2D(verts_2D)            
            all_points = np.concatenate((verts_2D,offset),axis=0)
            try:
                total_hull = ConvexHull(all_points,qhull_options='Pp') #ignores very small angles
                total_area = self._PolyArea2D(all_points[total_hull.vertices])
            except sp.spatial.qhull.QhullError:
                print(all_points)
                total_area =999

            if (total_area>original_area): # Throat is fully occluded
                " Don't do anything "
                Error = 5
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
    
    def _outer_offset(self,verts,fibre_rad):
        r"""
        Routine to loop through all verts and calculate offset position based on neighbours either side. Verts must be in hull order
        """    
        offset = []
        for i,vert in enumerate(verts):
            " Collect three adjacent points and compute the offset of the first "
            triplet = (vert, np.roll(verts,-1,axis=0)[i],np.roll(verts,1,axis=0)[i])
            offset.append(self._offset_vertex(triplet,fibre_rad))
        offset = np.asarray(offset)
        
        return offset    
    
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
    
        if alpha == 0: # this may cause problems in terms of which way to offset!!!
            #x = rad*np.cos(theta)
            #y = rad*np.sin(theta)
            x=0
            y=0
        else:
            x = rad*np.cos(alpha+theta)/np.sin(alpha)
            y = rad*np.sin(alpha+theta)/np.sin(alpha)

        "Add the midpoint back in"
        output = [x+p0[0],y+p0[1]]    
    
        return output

    def _dist2(self,p1, p2):
        r"""
        Pythagoras to compute the square of the distance between two points (in 2D)
        """
        return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2

    def _fuse(self,points, d):
        r"""
        Fuse points together wihin a certain range
        !!!! Not used anymore - superceded by line methods !!!!
        """
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
        r"""
        Work out the span of the points and therefore the range for fusing them together then call fuse
        !!!! Not used anymore - superceded by line methods !!!!
        """
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
        r"""
        returns the area of a 2D polygon given the set of points defining the convex hull in correct order
        """
        lines = np.hstack([pts,np.roll(pts,-1,axis=0)])
        area = 0.5*abs(sum(x1*y2-x2*y1 for x1,y1,x2,y2 in lines))
        return area

    def _PolyPerimeter2D(self,pts):
        r"""
        returns the perimeter of a 2D polygon given the set of points defining the convex hull in correct order
        """
        lines = np.hstack([pts,np.roll(pts,-1,axis=0)])
        perimeter = sum(np.sqrt((x2-x1)**2+(y2-y1)**2) for x1,y1,x2,y2 in lines)
        return perimeter
    
    def _get_hull_volume(self,points):
        
        r"""
        Calculate the volume of a set of points by dividing the bounding surface into triangles and working out the volume of all the pyramid elements
        connected to the volume centroid
        """
        " remove any duplicate points - this messes up the triangulation "        
        points = np.asarray(self._unique_list(points))       
        try:
            tri = Delaunay(points)
        except sp.spatial.qhull.QhullError:
            print(points)
        " We only want points included in the convex hull to calculate the centroid "
        #hull_points = np.unique(tri.convex_hull)#could technically use network pore centroids here but this function may be called at other times
        hull_centroid = sp.array([points[:,0].mean(),points[:,1].mean(),points[:,2].mean()])
        hull_volume = 0.0
        for ia, ib, ic in tri.convex_hull:
            " Points making each triangular face "
            " Collection of co-ordinates of each point in this face "
            face_x = points[[ia,ib,ic]][:,0]
            face_y = points[[ia,ib,ic]][:,1]
            face_z = points[[ia,ib,ic]][:,2]
            " Average of each co-ordinate is the centroid of the face "
            face_centroid = [face_x.mean(),face_y.mean(),face_z.mean()]
            face_centroid_vector = face_centroid - hull_centroid
            " Vectors of the sides of the face used to find normal vector and area "
            vab = points[ib] - points[ia]
            vac = points[ic] - points[ia]
            vbc = points[ic] - points[ib] # used later for area
            " As vectors are co-planar the cross-product will produce the normal vector of the face "
            face_normal = np.cross(vab,vac)
            face_unit_normal = face_normal/np.linalg.norm(face_normal)
            " As triangles are orientated randomly in 3D we could either transform co-ordinates to align with a plane and perform 2D operations "
            " to work out the area or we could work out the lengths of each side and use Heron's formula which is easier"
            " Using Delaunay traingulation will always produce triangular faces but if dealing with other polygons co-ordinate transfer may be necessary "
            a = np.linalg.norm(vab)
            b = np.linalg.norm(vbc)
            c = np.linalg.norm(vac)
            " Semiperimeter "
            s = 0.5*(a+b+c)
            face_area = np.sqrt(s*(s-a)*(s-b)*(s-c))
            " Now the volume of the pyramid section defined by the 3 face points and the hull centroid can be calculated "
            pyramid_volume = np.abs(np.dot(face_centroid_vector,face_unit_normal)*face_area/3)
            " Each pyramid is summed together to calculate the total volume "
            hull_volume += pyramid_volume
    
        return hull_volume
    
    def _unique_list(self,input_list):
        r"""
        For a given list (of points) remove any duplicates
        """
        output_list = []
        dim = np.shape(input_list)[1]
        for i in input_list:
            match=False
            for j in output_list:
                if dim == 3:
                    if (i[0]==j[0]) and (i[1]==j[1]) and (i[2]==j[2]):
                        match=True
                elif dim == 2:
                    if (i[0]==j[0]) and (i[1]==j[1]):
                        match=True
                elif dim ==1:
                    if (i[0]==j[0]):
                        match=True
            if match==False:
                output_list.append(i)
        return output_list
        
    def _symmetric_difference(self,list_a,list_b):
        r"""
        Return the combination of two lists without common elements (necessary as sets cannot contain mutable objects)
        """
        sym_diff=[]
        sorted_list_a = np.sort(list_a)
        sorted_list_b = np.sort(list_b)
        " Add elements in list a if not in list b "
        for element_a in sorted_list_a:
            match = False
            for element_b in sorted_list_b:
                if all(element_a==element_b):
                    match = True
            if match==False:
                sym_diff.append(element_a)
        " Add elements in list b if not in list a "
        for element_b in sorted_list_b:
            match = False
            for element_a in sorted_list_a:
                if all(element_a==element_b):
                    match = True
            if match==False:
                sym_diff.append(element_b)
        return sym_diff
        
    def _new_point(self,pairs):
        " Passed 2 pairs of points defining lines either side of overlapped offset vertices "
        " need to calculate the new point to offset from given the orientation of the outer fibres "

        m1,c1 = self._line_equation(pairs[0])
        m2,c2 = self._line_equation(pairs[1])

        if (m1 == np.inf):
            "line1 is a straight line x=c1"
            x=c1
            y=(m2*c1)+c2
        elif (m2 == np.inf):
            "line2 is a straight line x=c2"
            x=c2
            y=(m1*c2)+c1
        else:
            x=(c2-c1)/(m1-m2)
            y=(m1*c2 - m2*c1)/(m1-m2)

        return x,y
    
    def _line_equation(self,points):
        r"""
        Return the gradient and x intercept of a straight line given 2 points
        """
        x_coords, y_coords = zip(*points)
        dy = y_coords[1]-y_coords[0]
        dx = x_coords[1]-x_coords[0]
        if dx==0:
            m=np.inf
            c=x_coords[1]
        else:
            m=dy/dx
            c=y_coords[1] - m*x_coords[1]

        return m,c
        
    def _set_overlap(self,verts,offset):
        r"""
        Given a set of vertices and a set of offset vertices, evaluate whether any of the offset vertices overlap
        This is then used to recalculate points from which to offset
        """    
        dim = len(verts)
        overlap_array = np.zeros(shape=(dim,dim))
        sweep_radius = np.zeros(len(verts))
        for i,zone_centre in enumerate(verts):
            sweep_radius[i] = np.sqrt(self._dist2(zone_centre,offset[i]))
            for j,test_point in enumerate(offset):
                test_radius = np.sqrt(self._dist2(zone_centre,test_point))
                if (test_radius < sweep_radius[i]):
                    overlap_array[i][j]=1
                    overlap_array[j][i]=1 # Fill in both so that they are both recalculated later i overlapping j doesn't necessarily mean j overlaps i
    
        line_points=self._line_points(overlap_array)
    
        return overlap_array,sweep_radius,line_points

    def _line_points(self,array):
        r"""
        We are passed a square array containing a list of overlap results. rows represent vertices and columns represent offset vertices
        If an overlap occurs in the span between offset j and offset i then a 1 will result in [i][j]
        As the vertices are in hull order and our aim is to create a new point from which to offset using connected fibres we want to
        identify the correct points to use to define our lines
        
          e__________ d  
           |        | c  
           |       /     
           |     /       
           |   /         
           |_/           
           a b           
        if we have 5 points in a hull a,b,c,d,e where a overlaps with b and c overlaps with d (and visa versa) the array will look like
        [0,1,0,0,0]
        [1,0,0,0,0]
        [0,0,0,1,0]
        [0,0,1,0,0]
        [0,0,0,0,0]
        
        This means that c and d are within a fibre's width of each other and the line between them does not represent the fibre
        Instead we want to extend the outer lines (bc and de) to see where they would meet and offset from this point.
        Roll up and down to find the first unoverlapped index from which to start each line from then go back one to get the two
        points to form a line.
        """
        dim=np.shape(array)[0]
        index=range(dim)
        line_points=np.ndarray(shape=[dim,dim],dtype=object)
        for i in range(dim):
            for j in range(dim):
                if array[i][j]==1 and array[j][i]==1:
                    " Roll forwards to find the first unoverlapped index"
                    k=1                    
                    while k < dim:
                        if np.roll(array[i],-k,axis=0)[j]==0:
                            break
                        else:
                            k+=1
                    " Save the indices of the first unoverlapped index and the one before to create the line "
                    forward_line = [np.roll(index,-k,axis=0)[j],np.roll(index,-(k-1),axis=0)[j]]
                    forward_line.sort()
                    " Roll backwards to find the first unoverlapped index "
                    k=1
                    while k < dim:
                        if np.roll(array[i],k,axis=0)[j]==0:
                            break
                        else:
                            k+=1
                    " Save the indices of the first unoverlapped index and the one before to create the line "
                    backward_line = [np.roll(index,k,axis=0)[j],np.roll(index,(k-1),axis=0)[j]]
                    backward_line.sort()
                    line_points[i][j]=(forward_line,backward_line)
    
        return line_points  
        
    def _all_overlap(self,array):
        r""" 
        Find out whether all offset vertices (columns) are overlapped by at least one other
        If so then throat is fully occluded 
        """
        dim=np.shape(array)[0]
        overlap=[False]*dim
        all_overlap=False
        for i in range(dim):
            for j in range(dim):
                if array[j][i]==1:
                    overlap[i]=True
        if sum(overlap)==dim:
            all_overlap = True
        
        return all_overlap
        
if __name__ == '__main__':
    pn = OpenPNM.Network.Delaunay(name='test_net')
    pn.generate(num_pores=100, domain_size=[0.0001,0.0001,0.0001],add_boundaries=True)
    test = OpenPNM.Geometry.Voronoi(loglevel=10,name='test_geom',locations=[0],network=pn)
    test.set_locations(pores=pn.pores('internal'),throats='all') # Don't really know what this does but is needed
    pn.regenerate_geometries()