r"""
===============================================================================
Submodule -- pore_volume
===============================================================================

"""
import scipy as sp

def _one_unique_list(input_list):
    output_list = []
    for i in input_list:
        match=False
        for j in output_list:
            if (i[0]==j[0]) and (i[1]==j[1]) and (i[2]==j[2]):
                match=True
        if match==False:
            output_list.append(i)
    return output_list

def _centroid(points):
    import numpy as np
    columns = np.shape(points)[1]
    if columns == 2:
        cen = np.array([points[:,0].mean(),points[:,1].mean()])
    elif columns == 3:
        cen = np.array([points[:,0].mean(),points[:,1].mean(),points[:,2].mean()])
    else:
        cen = 0
    return cen

def _get_hull_volume(points):
    import numpy as np
    from scipy.spatial import Delaunay

    #points = np.random.randn(20, 3) # Random points in space
    #points = np.array([[0,0,0],[0,1,0],[1,0,0],[1,1,0],[0.5,0.5,3]]) # Regular Conical Pyramid - Volume = 1.0
    try:
        tri = Delaunay(points)
    except sp.spatial.qhull.QhullError:
        print(points)
    " We only want points included in the convex hull to calculate the centroid "
    hull_points = np.unique(tri.convex_hull)
    #hull_centroid = np.array([points[hull_points,0].mean(),points[hull_points,1].mean(),points[hull_points,2].mean()])
    hull_centroid = _centroid(points[hull_points])
    # -- Make a list of faces, [(p1, p2, p3), ...];  pj = (xj, yj, zj)

    #faces = [] # list of points making each face
    #face_centers = [] # co-ordinates of the face centroids
    #face_normals = [] # normal vector of the faces
    #face_areas = [] # Area of each face
    hull_volume = 0.0
    for ia, ib, ic in tri.convex_hull:
        " Points making each triangular face "
        #faces.append(points[[ia, ib, ic]])
        " Collection of co-ordinates of each point in this face "
        face_x = points[[ia,ib,ic]][:,0]
        face_y = points[[ia,ib,ic]][:,1]
        face_z = points[[ia,ib,ic]][:,2]
        " Average of each co-ordinate is the centroid of the face "
        face_centroid = [face_x.mean(),face_y.mean(),face_z.mean()]
        #face_centers.append(face_centroid)
        face_centroid_vector = face_centroid - hull_centroid
        " Vectors of the sides of the face used to find normal vector and area "
        vab = points[ib] - points[ia]
        vac = points[ic] - points[ia]
        vbc = points[ic] - points[ib] # used later for area
        " As vectors are co-planar the cross-product will produce the normal vector of the face "
        face_normal = np.cross(vab,vac)
        face_unit_normal = face_normal/np.linalg.norm(face_normal)
        #face_normals.append(face_unit_normal)
        " As triangles are orientated randomly in 3D we could either transform co-ordinates to align with a plane and perform 2D operations "
        " to work out the area or we could work out the lengths of each side and use Heron's formula which is easier"
        " Using Delaunay traingulation will always produce triangular faces but if dealing with other polygons co-ordinate transfer may be necessary "
        a = np.linalg.norm(vab)
        b = np.linalg.norm(vbc)
        c = np.linalg.norm(vac)
        " Semiperimeter "
        s = 0.5*(a+b+c)
        face_area = np.sqrt(s*(s-a)*(s-b)*(s-c))
        #face_areas.append(face_area)
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
    Calculate volume from the convex hull of the vertices making the Voronoi Polyhedron
    """
    all_verts = network.get_pore_data(prop='vertices')
    value = sp.ndarray(len(all_verts))
    for i,pore_verts in enumerate(all_verts):
        if pore_verts != "unbounded":
            value[i] = _get_hull_volume(pore_verts)
        else:
            value[i] = 0.0
    network.set_data(prop=propname,pores=geometry.pores(),data=value)
    
def voronoi2(geometry,
            network,
            propname,
            **params):
    r"""
    Calculate volume from the convex hull of the offset vertices making the throats
    """
    conns = network.get_throat_data(prop='connections')
    verts = network.get_throat_data(prop='offset_verts') #This won't work as pore data is done before throat, need to calculate offset verts at network generation stage 
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
        if num_connections > 1:#need at least 4 sides to create a bound volume (need to think about how occluded faces are treated here)
            unique_throat_list=_one_unique_list(throat_vert_list)
            throat_array=sp.asarray(unique_throat_list)
            value[my_pore]=_get_hull_volume(throat_array)
        else:
            value[my_pore]=0.0

    network.set_data(prop=propname,pores=geometry.pores(),data=value)