r"""
===============================================================================
Submodule -- throat_area
===============================================================================

"""
import scipy as sp
import scipy.stats as spst
import numpy as np
import _transformations as tr
from scipy.spatial import ConvexHull
from math import atan2

def offset_vertex(points,rad = 0.01):
     
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

def dist2(p1, p2):
    return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2

def fuse(points, d):
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
                if dist2(points[i], points[j]) < d2:
                    point[0] += points[j][0]
                    point[1] += points[j][1]
                    count+=1
                    taken[j] = True
            point[0] /= count
            point[1] /= count
            ret.append((point[0], point[1]))
    return ret

def fuse_verts(verts,percentage=0.05):
    #Work out largest span
    x_span = max(verts[:,0])- min(verts[:,0])
    y_span = max(verts[:,1])- min(verts[:,1])
    if x_span > y_span:
        tolerance = x_span*percentage
    else:
        tolerance = y_span*percentage
    #fuse vertices lying within 5% of the largest span    
    return fuse(verts,tolerance)

def PolyArea2D(pts):
    lines = np.hstack([pts,np.roll(pts,-1,axis=0)])
    area = 0.5*abs(sum(x1*y2-x2*y1 for x1,y1,x2,y2 in lines))
    return area

def PolyPerimeter2D(pts):
    lines = np.hstack([pts,np.roll(pts,-1,axis=0)])
    perimeter = sum(np.sqrt((x2-x1)**2+(y2-y1)**2) for x1,y1,x2,y2 in lines)
    return perimeter
    
def _get_throat_geom(verts,normal):
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
    fused_verts=fuse_verts(verts_2D)
    if len(fused_verts) <3:
        #we fused too many
        fused_verts=fuse_verts(verts_2D,0.0025)
    offset = []
    fibre_rad = 1e-06
    for i,vert in enumerate(fused_verts):
        " Collect three adjacent points and compute the offset of the first "
        triplet = (vert, np.roll(fused_verts,-1,axis=0)[i],np.roll(fused_verts,1,axis=0)[i])
        offset.append(offset_vertex(triplet,fibre_rad))
    offset = np.asarray(offset)
    
    original_area = PolyArea2D(verts_2D)            
    all_points = np.concatenate((verts_2D,offset),axis=0)
    try:
        total_hull = ConvexHull(all_points,qhull_options='Pp') #ignores very small angles
        total_area = PolyArea2D(all_points[total_hull.vertices])
    except sp.spatial.qhull.QhullError:
        print(all_points)
        total_area =999
    #total_area = PolyArea2D(all_points[total_hull.vertices])
    if (total_area>original_area): # Throat is fully occluded
        throat_area = 0.0
        throat_perimeter = 0.0
        output_offset = []
    else:
        offset_hull = ConvexHull(offset)
        offset_verts_2D = offset[offset_hull.vertices]
        throat_area = PolyArea2D(offset_verts_2D)
        throat_perimeter = PolyPerimeter2D(offset_verts_2D)
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
    
def constant(geometry,
             network,
             propname,
             value,
             **params):
    r"""
    Assigns specified constant value
    """
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

def cylinder(geometry,
             network,
             propname,
             diameter='diameter',
             **params):
    r"""
    Calculate throat area for a cylindrical throat
    """
    D = network.get_data(prop=diameter,throats=geometry.throats())
    value = sp.constants.pi/4*(D)**2
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

def cuboid(geometry,
           network,
           propname,
           diameter='diameter',
           **params):
    r"""
    Calculate throat area for a cuboid throat
    """
    D = network.get_data(prop=diameter,throats=geometry.throats())
    value = (D)**2
    network.set_data(prop=propname,throats=geometry.throats(),data=value)
    
def voronoi(geometry,
            network,
            propname,
            **params):
    r"""
    Calculate the area of the voronoi facet shared between the connecting pores
    - Later add the vertex offset routine to adjust for the fibre occlusion
    - This step must be performed first and then subsequently diameter could be calculated
    """
    connections = network.get_throat_data(prop='connections')
    coords = network.get_pore_data(prop='coords')
    verts = network.get_pore_data(prop='vertices')
    normals = coords[connections[:,0]]-coords[connections[:,1]]
    area = sp.ndarray(len(connections),dtype=object)
    perimeter = sp.ndarray(len(connections),dtype=object)
    centroids = sp.ndarray(len(connections),dtype=object)
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
            area[i],perimeter[i],offset_verts[i] = _get_throat_geom(shared_verts,normals[i])
            centroids[i]=network._centroid(shared_verts)
        #else:
            #area[i] = 0.0
            #perimeter = 0.0
            #centroids[i]=0.0
    network.set_data(prop=propname,throats=geometry.throats(),data=area)
    network.set_data(prop='perimeter',throats=geometry.throats(),data=perimeter)
    network.set_data(prop='centroid',throats=geometry.throats(),data=centroids)
    network.set_data(prop='offset_verts',throats=geometry.throats(),data=offset_verts)
    " Temporary Code to remove the smallest 2% of throat area connections "
    " This can be used in algorithms to ignore certain connections if required - like in the range of capillary pressures in OP "
    smallest_throat = min(area)
    largest_throat = max(area)
    remove_range = smallest_throat + ((largest_throat-smallest_throat)*0.02)
    excluded_throats = []
    for i,throat in enumerate(connections):
        if area[i]<=remove_range:
        #if area[i]==0.0:
            excluded_throats.append(i)
    excluded_throats = np.asarray(excluded_throats)
    if len(excluded_throats) > 0:
        network.trim(throats=excluded_throats)

    