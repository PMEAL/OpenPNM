r"""
===============================================================================
pore_volume --
===============================================================================

"""
import scipy as _sp
import numpy as np
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
import OpenPNM.Utilities.misc as misc
from scipy import ndimage
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


def inhull(geometry, xyz, pore, tol=1e-7):
    r"""
    Tests whether points lie within a convex hull or not.
    Computes a tesselation of the hull works out the normals of the facets.
    Then tests whether dot(x.normals) < dot(a.normals) where a is the the
    first vertex of the facets
    """
    xyz = np.around(xyz, 10)
    # Work out range to span over for pore hull
    xmin = xyz[:, 0].min()
    xr = (np.ceil(xyz[:, 0].max())-np.floor(xmin)).astype(int)+1
    ymin = xyz[:, 1].min()
    yr = (np.ceil(xyz[:, 1].max())-np.floor(ymin)).astype(int)+1
    zmin = xyz[:, 2].min()
    zr = (np.ceil(xyz[:, 2].max())-np.floor(zmin)).astype(int)+1

    origin = np.array([xmin, ymin, zmin])
    # start index
    si = np.floor(origin).astype(int)
    xyz -= origin
    dom = np.zeros([xr, yr, zr], dtype=np.uint8)
    indx, indy, indz = np.indices((xr, yr, zr))
    # Calculate the tesselation of the points
    hull = ConvexHull(xyz)
    # Assume 3d for now
    # Calculate normals from the vector cross product of the vectors defined
    # by joining points in the simplices
    vab = xyz[hull.simplices[:, 0]]-xyz[hull.simplices[:, 1]]
    vac = xyz[hull.simplices[:, 0]]-xyz[hull.simplices[:, 2]]
    nrmls = np.cross(vab, vac)
    # Scale normal vectors to unit length
    nrmlen = np.sum(nrmls**2, axis=-1)**(1./2)
    nrmls = nrmls*np.tile((1/nrmlen), (3, 1)).T
    # Center of Mass
    center = np.mean(xyz, axis=0)
    # Any point from each simplex
    a = xyz[hull.simplices[:, 0]]
    # Make sure all normals point inwards
    dp = np.sum((np.tile(center, (len(a), 1))-a)*nrmls, axis=-1)
    k = dp < 0
    nrmls[k] = -nrmls[k]
    # Now we want to test whether dot(x,N) >= dot(a,N)
    aN = np.sum(nrmls*a, axis=-1)

    for plane_index in range(len(a)):
        eqx = nrmls[plane_index][0]*(indx)
        eqy = nrmls[plane_index][1]*(indy)
        eqz = nrmls[plane_index][2]*(indz)
        xN = eqx + eqy + eqz
        dom[xN - aN[plane_index] >= 0-tol] += 1
    dom[dom < len(a)] = 0
    dom[dom == len(a)] = 1
    ds = np.shape(dom)
    temp_arr = np.zeros_like(geometry._hull_image, dtype=bool)
    temp_arr[si[0]:si[0]+ds[0], si[1]:si[1]+ds[1], si[2]:si[2]+ds[2]] = dom
    geometry._hull_image[temp_arr] = pore
    hull_num = np.sum(dom)
    dom = dom * geometry._fibre_image[si[0]:si[0]+ds[0], si[1]:si[1]+ds[1],
                                      si[2]:si[2]+ds[2]]
    pore_num = np.sum(dom)
    fibre_num = hull_num - pore_num

    del temp_arr
    return pore_num, fibre_num


def _voxel_centroid(image, pores=None, vox_len=1e-6):
    r"""
    Calculate Pore Centroid from indices of the voronoi voxel image generated
    in _get_voxel_volume
    """
    if pores is None:
        pores = np.unique(image)
    centroids = np.zeros([len(pores), 3])
    for i, pore in enumerate(pores):
        px, py, pz = np.where(image == pore)
        try:
            # If one is empty then all will be - no pore volume for this pore
            if len(px) > 0:
                cx = np.mean(px)
                cy = np.mean(py)
                cz = np.mean(pz)
            centroids[i] = np.array([cx, cy, cz])*vox_len
        except:
            logger.warn("Some centroid data may be invalid, look for zeros")
    return centroids


def _get_vertex_range(verts):
    # Find the extent of the vetrices
    [vxmin, vxmax, vymin, vymax, vzmin, vzmax] = [1e20, 0, 1e20, 0, 1e20, 0]
    for vert in verts:
        if np.min(vert[:, 0]) < vxmin:
            vxmin = np.min(vert[:, 0])
        if np.max(vert[:, 0]) > vxmax:
            vxmax = np.max(vert[:, 0])
        if np.min(vert[:, 1]) < vymin:
            vymin = np.min(vert[:, 1])
        if np.max(vert[:, 1]) > vymax:
            vymax = np.max(vert[:, 1])
        if np.min(vert[:, 2]) < vzmin:
            vzmin = np.min(vert[:, 2])
        if np.max(vert[:, 2]) > vzmax:
            vzmax = np.max(vert[:, 2])
    return [vxmin, vxmax, vymin, vymax, vzmin, vzmax]


def _get_fibre_image(network, cpores, vox_len, fibre_rad):
    r"""
    Produce image by filling in voxels along throat edges using Bresenham line
    Then performing distance transform on fibre voxels to erode the pore space
    """

    cthroats = network.find_neighbor_throats(pores=cpores)

    # Below method copied from geometry model throat.vertices
    # Needed now as network may not have all throats assigned to geometry
    # i.e network['throat.vertices'] could return garbage
    verts = _sp.ndarray(network.num_throats(), dtype=object)
    for i in range(len(verts)):
        verts[i] = _sp.asarray(list(network["throat.vert_index"][i].values()))
    cverts = verts[cthroats]
    [vxmin, vxmax, vymin, vymax, vzmin, vzmax] = _get_vertex_range(cverts)
    # Translate vertices so that minimum occurs at the origin
    for index in range(len(cverts)):
        cverts[index] -= np.array([vxmin, vymin, vzmin])
    # Find new size of image array
    cdomain = np.around(np.array([(vxmax-vxmin),
                                  (vymax-vymin),
                                  (vzmax-vzmin)]), 6)
    logger.info("Creating fibre domain range: " + str(np.around(cdomain, 5)))
    lx = np.int(np.around(cdomain[0]/vox_len)+1)
    ly = np.int(np.around(cdomain[1]/vox_len)+1)
    lz = np.int(np.around(cdomain[2]/vox_len)+1)
    # Try to create all the arrays we will need at total domain size
    try:
        pore_space = np.ones([lx, ly, lz], dtype=np.uint8)
        fibre_space = np.zeros(shape=[lx, ly, lz], dtype=np.uint8)
        dt = np.zeros([lx, ly, lz], dtype=float)
        # Only need one chunk
        cx = cy = cz = 1
        chunk_len = np.max(np.shape(pore_space))
    except:
        logger.info("Domain too large to fit into memory so chunking domain"
                    "to process image, this may take some time")
        # Do chunking
        chunk_len = 100
        if (lx > chunk_len):
            cx = np.ceil(lx/chunk_len).astype(int)
        else:
            cx = 1
        if (ly > chunk_len):
            cy = np.ceil(ly/chunk_len).astype(int)
        else:
            cy = 1
        if (lz > chunk_len):
            cz = np.ceil(lz/chunk_len).astype(int)
        else:
            cz = 1

    # Get image of the fibres
    line_points = bresenham(cverts, vox_len/2)
    line_ints = (np.around((line_points/vox_len), 0)).astype(int)
    for x, y, z in line_ints:
        try:
            pore_space[x][y][z] = 0
        except IndexError:
            logger.warning("Some elements in image processing are out" +
                           "of bounds")

    num_chunks = np.int(cx*cy*cz)
    cnum = 1
    for ci in range(cx):
        for cj in range(cy):
            for ck in range(cz):
                # Work out chunk range
                logger.info("Processing Fibre Chunk: "+str(cnum)+" of " +
                            str(num_chunks))
                cxmin = ci*chunk_len
                cxmax = np.int(np.ceil((ci+1)*chunk_len + 5*fibre_rad))
                cymin = cj*chunk_len
                cymax = np.int(np.ceil((cj+1)*chunk_len + 5*fibre_rad))
                czmin = ck*chunk_len
                czmax = np.int(np.ceil((ck+1)*chunk_len + 5*fibre_rad))
                # Don't overshoot
                if cxmax > lx:
                    cxmax = lx
                if cymax > ly:
                    cymax = ly
                if czmax > lz:
                    czmax = lz
                dt = ndimage.distance_transform_edt(pore_space[cxmin:cxmax,
                                                               cymin:cymax,
                                                               czmin:czmax])
                fibre_space[cxmin:cxmax,
                            cymin:cymax,
                            czmin:czmax][dt <= fibre_rad] = 0
                fibre_space[cxmin:cxmax,
                            cymin:cymax,
                            czmin:czmax][dt > fibre_rad] = 1
                del dt
                cnum += 1
    del pore_space
    return fibre_space


def bresenham(faces, dx):
    line_points = []
    for face in faces:
        # Get in hull order
        fx = face[:, 0]
        fy = face[:, 1]
        fz = face[:, 2]
        # Find the axis with the smallest spread and remove it to make 2D
        if (np.std(fx) < np.std(fy)) and (np.std(fx) < np.std(fz)):
            f2d = np.vstack((fy, fz)).T
        elif (np.std(fy) < np.std(fx)) and (np.std(fy) < np.std(fz)):
            f2d = np.vstack((fx, fz)).T
        else:
            f2d = np.vstack((fx, fy)).T
        hull = ConvexHull(f2d, qhull_options='QJ Pp')
        face = np.around(face[hull.vertices], 6)
        for i in range(len(face)):
            vec = face[i]-face[i-1]
            vec_length = np.linalg.norm(vec)
            increments = np.ceil(vec_length/dx)
            check_p_old = np.array([-1, -1, -1])
            for x in np.linspace(0, 1, increments):
                check_p_new = face[i-1]+(vec*x)
                if np.sum(check_p_new - check_p_old) != 0:
                    line_points.append(check_p_new)
                    check_p_old = check_p_new
    return np.asarray(line_points)


def sphere(geometry, pore_diameter='pore.diameter', **kwargs):
    r"""
    Calculate pore volume from diameter assuming a spherical pore body

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key to the pore diameter values

    """
    diams = geometry[pore_diameter]
    value = _sp.pi/6*diams**3
    return value


def cube(geometry, pore_diameter='pore.diameter', **kwargs):
    r"""
    Calculate pore volume from diameter assuming a cubic pore body

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key to the pore diameter values

    """
    diams = geometry[pore_diameter]
    value = diams**3
    return value


def _get_hull_volume(points):
    r"""
    Calculate the volume of a set of points by dividing the bounding surface
    into triangles and working out the volume of all the pyramid elements
    connected to the volume centroid
    """
    # Remove any duplicate points - this messes up the triangulation
    points = _sp.asarray(misc.unique_list(np.around(points, 10)))
    try:
        tri = Delaunay(points, qhull_options='QJ Pp')
    except _sp.spatial.qhull.QhullError:
        logger.error("Volume suspect for points: " + str(points))
    # We only want points included in the convex hull to calculate the centroid
    hull_centroid = _sp.array([points[:, 0].mean(),
                               points[:, 1].mean(),
                               points[:, 2].mean()])
    hull_volume = 0.0
    pyramid_COMs = []
    for ia, ib, ic in tri.convex_hull:
        # Points making each triangular face
        # Collection of co-ordinates of each point in this face
        face_x = points[[ia, ib, ic]][:, 0]
        face_y = points[[ia, ib, ic]][:, 1]
        face_z = points[[ia, ib, ic]][:, 2]
        # Average of each co-ordinate is the centroid of the face
        face_centroid = [face_x.mean(), face_y.mean(), face_z.mean()]
        face_centroid_vector = face_centroid - hull_centroid
        # Vectors of the sides of the face used to find normal vector and area
        vab = points[ib] - points[ia]
        vac = points[ic] - points[ia]
        vbc = points[ic] - points[ib]
        # As vectors are co-planar the cross-product will produce the normal
        # vector of the face
        face_normal = _sp.cross(vab, vac)
        try:
            face_unit_normal = face_normal/_sp.linalg.norm(face_normal)
        except RuntimeWarning:
            logger.error('Pore Volume Error:' + str(vab) + ' ' + str(vac))
        """
        As triangles are orientated randomly in 3D we could either transform
        co-ordinates to align with a plane and perform 2D operations to work
        out the area or we could work out the lengths of each side and use
        Heron's formula which is easier. Using Delaunay traingulation will
        always produce triangular faces but if dealing with other polygons
        co-ordinate transfer may be necessary
        """
        a = _sp.linalg.norm(vab)
        b = _sp.linalg.norm(vbc)
        c = _sp.linalg.norm(vac)
        # Semiperimeter
        s = 0.5*(a + b + c)
        face_area = _sp.sqrt(s*(s-a)*(s-b)*(s-c))
        # Now the volume of the pyramid section defined by the 3 face points
        # and the hull centroid can be calculated "
        pyramid_volume = _sp.absolute(_sp.dot(face_centroid_vector,
                                              face_unit_normal)*face_area/3)
        # Each pyramid is summed together to calculate the total volume
        hull_volume += pyramid_volume
        # The Centre of Mass will not be the same as the geometrical centroid.
        # Weighted adjustment is calculated from pyramid centroid and volume.
        vha = points[ia] - hull_centroid
        vhb = points[ib] - hull_centroid
        vhc = points[ic] - hull_centroid
        pCOM = ((vha+vhb+vhc)/4)*pyramid_volume
        pyramid_COMs.append(pCOM)
    if _sp.isnan(hull_volume):
        hull_volume = 0.0
    if hull_volume > 0:
        hull_COM = hull_centroid + _sp.mean(_sp.asarray(pyramid_COMs),
                                            axis=0) / hull_volume
    else:
        hull_COM = hull_centroid

    return hull_volume, hull_COM


def voronoi(network, geometry, **kwargs):
    r"""
    Calculate volume from the convex hull of the offset vertices making the
    throats surrounding the pore. Also calculate the centre of mass for the
    volume.
    """
    pores = geometry.map_pores(network, geometry.pores())
    Np = len(pores)
    volume = _sp.zeros(Np)
    com = _sp.zeros([Np, 3])
    for i in range(Np):
        throat_vert_list = []
        net_throats = network.find_neighbor_throats([pores[i]])
        geom_throats = network.map_throats(target=geometry,
                                           throats=net_throats,
                                           return_mapping=True)['target']
        if len(geom_throats) > 1:
            for throat in geom_throats:
                geom_throat_verts = geometry["throat.offset_vertices"][throat]
                if geom_throat_verts is not None:
                    for j in range(len(geom_throat_verts)):
                        throat_vert_list.append(geom_throat_verts[j])
            throat_array = _sp.asarray(throat_vert_list)
            if len(throat_array) > 4:
                volume[i], com[i] = _get_hull_volume(throat_array)
            else:
                volume[i] = 0
        elif len(geom_throats) == 1 and 'throat.centroid' in geometry.props():
                com[i] = geometry['throat.centroid'][geom_throats]
                volume[i] = 0
    # Find any pores with centroids at origin and use the mean of the pore
    # vertices.  Not doing this messes up hydraulic conductances using centre
    # to centre
    ps = np.where(~com.any(axis=1))[0]
    if len(ps) > 0:
        for pore in ps:
            com[pore] = np.mean(geometry['pore.vertices'][pore], axis=0)
    geometry['pore.centroid'] = com

    return volume


def in_hull_volume(network, geometry, fibre_rad, vox_len=1e-6, **kwargs):
    r"""
    Work out the voxels inside the convex hull of the voronoi vertices of each
    pore
    """
    Np = network.num_pores()
    geom_pores = geometry.map_pores(network, geometry.pores())
    volume = _sp.zeros(Np)
    pore_vox = _sp.zeros(Np, dtype=int)
    fibre_vox = _sp.zeros(Np, dtype=int)
    voxel = vox_len**3
    try:
        nbps = network.pores('boundary', mode='not')
    except KeyError:
        # Boundaries have not been generated
        nbps = network.pores()
    # Voxel length
    fibre_rad = np.around((fibre_rad-(vox_len/2))/vox_len, 0).astype(int)

    # Get the fibre image
    fibre_image = _get_fibre_image(network, geom_pores, vox_len, fibre_rad)
    # Save as private variables
    geometry._fibre_image = fibre_image
    hull_image = np.ones_like(fibre_image, dtype=np.uint16)*-1
    geometry._hull_image = hull_image
    for pore in nbps:
        logger.info("Processing Pore: "+str(pore+1)+" of "+str(len(nbps)))
        verts = np.asarray([i for i in network["pore.vert_index"][pore].values()])
        verts = np.asarray(misc.unique_list(np.around(verts, 6)))
        verts /= vox_len
        pore_vox[pore], fibre_vox[pore] = inhull(geometry, verts, pore)

    volume = pore_vox*voxel
    geometry["pore.fibre_voxels"] = fibre_vox[geom_pores]
    geometry["pore.pore_voxels"] = pore_vox[geom_pores]

    return volume[geom_pores]
