r"""
===============================================================================
pore_centroid -- 
===============================================================================

"""
import scipy as _sp
import numpy as np
import OpenPNM.Utilities.misc as misc
from scipy.spatial import Delaunay

def _get_hull_com(points):
    r"""
    Divide the hull of the offset vertices into pyramids and work out their centre of masses then 
    calculate a weighted COM for the .
    N.B This is inefficient as all work done here to calculate COM could be used to calculate pore volume too 
        But processes may be separate and alternative volume or centroid methods maybe used instead
    """
    " remove any duplicate points - this messes up the triangulation "        
    points = _sp.asarray(misc.unique_list(points))       
    try:
        tri = Delaunay(points,qhull_options='QJ Pp')
    except _sp.spatial.qhull.QhullError:
        print(points)
    " We only want points included in the convex hull to calculate the centroid "
    hull_centroid = _sp.array([points[:,0].mean(),points[:,1].mean(),points[:,2].mean()])
    hull_volume = 0.0
    pyramid_COMs = []
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
        #face_COM = (vab+vac)/3
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
        " The Centre of Mass will not be the same as the geometrical centroid " 
        " A weighted adjustment can be calculated from the pyramid centroid and volume "
        vha = points[ia]-hull_centroid
        vhb = points[ib]-hull_centroid
        vhc = points[ic]-hull_centroid
        pCOM = ((vha+vhb+vhc)/4)*pyramid_volume
        pyramid_COMs.append(pCOM)
    
    if hull_volume>0:    
        hull_COM = hull_centroid + _sp.mean(_sp.asarray(pyramid_COMs),axis=0)/hull_volume
    else:
        hull_COM = hull_centroid
    
    return hull_COM

def voronoi(geometry,
            pore_vertices='pore.vertices',
            **kwargs):
    r"""
    Calculate the centroid of the pore from the voronoi vertices - C.O.M
    """
    #network = geometry._net    
    #pores = geometry['pore.map']
    verts = geometry[pore_vertices]
    value = _sp.ndarray([len(verts),3])
    for i,vert in enumerate(verts):
        value[i] = _sp.array([vert[:,0].mean(),vert[:,1].mean(),vert[:,2].mean()])
    return value

def voronoi2(geometry,
             vertices='throat.centroid',
             **kwargs):
    r"""
    Calculate the centroid from the mean of the throat centroids
    """
    value = _sp.ndarray([geometry.num_pores(),3])
    pore_map = geometry.map_pores(geometry.pores(),geometry._net)
    for geom_pore,net_pore in pore_map:
        net_throats = geometry._net.find_neighbor_throats(net_pore)
        geom_throats = geometry._net.map_throats(net_throats,geometry)[:,1]
        verts = geometry[vertices][geom_throats]
        value[geom_pore]=_sp.array([verts[:,0].mean(),verts[:,1].mean(),verts[:,2].mean()])
    return value

def centre_of_mass(network,
                   geometry,
                   **kwargs):
    r"""
    Calculate the centre of mass from the offset throat vertices
    """
    pores = geometry['pore.map']
    Np = len(pores)
    value = _sp.zeros([Np,3])
    for i in range(Np):
        throat_vert_list = []
        throats=network.find_neighbor_throats([pores[i]])
        if len(throats) > 1:        
            for throat in throats:
                try:
                    geom_throat = geometry['throat.map'].tolist().index(throat)
                    geom_throat_verts = geometry["throat.offset_vertices"][geom_throat]
                    
                    for j in range(len(geom_throat_verts)):
                        throat_vert_list.append(geom_throat_verts[j])
                except ValueError:
                    " Throat is not part of this geometry "
            throat_array=_sp.asarray(throat_vert_list)
            if len(throat_array)>4: #Why 4???
                value[i] = _get_hull_com(throat_array)
        elif len(throats) == 1 and 'throat.centroid' in geometry.props():
                geom_throat = geometry['throat.map'].tolist().index(throats)
                value[i]=geometry['throat.centroid'][geom_throat]
    "Find any pores with centroids at origin and use the mean of the pore vertices instead"
    "Not doing this messes up hydraulic conductances using centre to centre"
    ps = np.where(~value.any(axis=1))[0]
    if len(ps) >0:
        for pore in ps:
            value[pore]=np.mean(geometry["pore.vertices"][pore],axis=0)

    return value