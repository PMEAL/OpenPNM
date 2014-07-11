r"""
===============================================================================
Submodule -- pore_volume
===============================================================================

"""
import scipy as sp
import numpy as np
from scipy.spatial import Delaunay
import OpenPNM.Utilities.misc as misc

def _get_hull_volume(points):
        
    r"""
    Calculate the volume of a set of points by dividing the bounding surface into triangles and working out the volume of all the pyramid elements
    connected to the volume centroid
    """
    " remove any duplicate points - this messes up the triangulation "        
    points = np.asarray(misc.unique_list(points))       
    try:
        tri = Delaunay(points)
    except sp.spatial.qhull.QhullError:
        print(points)
    " We only want points included in the convex hull to calculate the centroid "
    #hull_points = np.unique(tri.convex_hull)#could technically use network pore centroids here but this function may be called at other times
    hull_centroid = np.array([points[:,0].mean(),points[:,1].mean(),points[:,2].mean()])
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
        face_normal = sp.cross(vab,vac)
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



def constant(geometry,
             network,
             propname,
             value,
             **params):
    r"""
    Assigns specified constant value
    """
    network.set_data(prop=propname,pores=geometry.pores(),data=value)

def sphere(geometry,
           network,
           propname,
           diameter='diameter',
           **params):
    r"""
    Calculate pore volume from diameter for a spherical pore body
    """
    value=sp.pi/6*network.get_data(prop=diameter,pores=geometry.pores())**3
    network.set_data(prop=propname,pores=geometry.pores(),data=value)
    
def cube(geometry,
         network,
         propname,
         diameter='diameter',
         **params):
    r"""
    Calculate pore volume from diameter for a cubic pore body
    """
    value=network.get_data(prop=diameter,pores=geometry.pores())**3
    network.set_data(prop=propname,pores=geometry.pores(),data=value)
    
def voronoi(geometry,
            network,
            propname,
            **params):
    r"""
    Calculate volume from the convex hull of the offset vertices making the throats
    """
    conns = network.get_throat_data(prop='conns')
    verts = network.get_throat_data(prop='offset_verts') 
    num_pores = network.num_pores()
    value = sp.ndarray(num_pores,dtype=object)
    for my_pore in range(num_pores):
        throat_vert_list = []
        num_connections = 0
        for idx,check_pores in enumerate(conns):
            if (check_pores[0] == my_pore) or (check_pores[1] == my_pore):
                num_connections +=1
                for vertex in range(len(verts[idx])):
                    throat_vert_list.append(verts[idx][vertex])
        if num_connections > 1:
            throat_array=sp.asarray(throat_vert_list)
            value[my_pore]= _get_hull_volume(throat_array)
        else:
            value[my_pore]=0.0

    network.set_data(prop=propname,pores=geometry.pores(),data=value)