r"""
===============================================================================
Submodule -- pore_volume
===============================================================================

"""
import scipy as _sp
from scipy.spatial import Delaunay
import OpenPNM.Utilities.misc as misc

def _get_hull_volume(points):
    r"""
    Calculate the volume of a set of points by dividing the bounding surface into triangles and working out the volume of all the pyramid elements
    connected to the volume centroid
    """
    " remove any duplicate points - this messes up the triangulation "        
    points = _sp.asarray(misc.unique_list(points))       
    try:
        tri = Delaunay(points)
    except _sp.spatial.qhull.QhullError:
        print(points)
    " We only want points included in the convex hull to calculate the centroid "
    hull_centroid = _sp.array([points[:,0].mean(),points[:,1].mean(),points[:,2].mean()])
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
        face_normal = _sp.cross(vab,vac)
        face_unit_normal = face_normal/_sp.linalg.norm(face_normal)
        " As triangles are orientated randomly in 3D we could either transform co-ordinates to align with a plane and perform 2D operations "
        " to work out the area or we could work out the lengths of each side and use Heron's formula which is easier"
        " Using Delaunay traingulation will always produce triangular faces but if dealing with other polygons co-ordinate transfer may be necessary "
        a = _sp.linalg.norm(vab)
        b = _sp.linalg.norm(vbc)
        c = _sp.linalg.norm(vac)
        " Semiperimeter "
        s = 0.5*(a+b+c)
        face_area = _sp.sqrt(s*(s-a)*(s-b)*(s-c))
        " Now the volume of the pyramid section defined by the 3 face points and the hull centroid can be calculated "
        pyramid_volume = _sp.absolute(_sp.dot(face_centroid_vector,face_unit_normal)*face_area/3)
        " Each pyramid is summed together to calculate the total volume "
        hull_volume += pyramid_volume
    
    return hull_volume

def sphere(geometry,
           pore_diameter='pore.diameter',
           **kwargs):
    r"""
    Calculate pore volume from diameter for a spherical pore body
    """
    diams = geometry[pore_diameter]
    value=_sp.pi/6*diams**3
    return value
    
def cube(geometry,
         pore_diameter='pore.diameter',
         **kwargs):
    r"""
    Calculate pore volume from diameter for a cubic pore body
    """
    diams = geometry[pore_diameter]
    value = diams**3
    return value
    
def voronoi(geometry,
            **kwargs):
    r"""
    Calculate volume from the convex hull of the offset vertices making the throats
    
    TODO: Optimise to only calculate volume of pores in geometry rather than selecting them afterwards
    """
    network = geometry._net    
    pores = geometry['pore.map']
    conns = network['throat.conns']
    verts = network['throat.offset_verts']
    Np = network.num_pores()
    value = _sp.ndarray(Np)
    for my_pore in range(Np):
        throat_vert_list = []
        num_connections = 0
        for idx,check_pores in enumerate(conns):
            if (check_pores[0] == my_pore) or (check_pores[1] == my_pore):
                num_connections +=1
                for vertex in range(len(verts[idx])):
                    throat_vert_list.append(verts[idx][vertex])
        if num_connections > 1:
            throat_array=_sp.asarray(throat_vert_list)
            value[my_pore]= _get_hull_volume(throat_array)
        else:
            value[my_pore]=0.0
    return value[pores]