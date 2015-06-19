r"""
===============================================================================
throat_offset_vertices -- Offeset throat vertices using a fibre radius parameter
===============================================================================

"""
import scipy as sp
import OpenPNM.Utilities.vertexops as vo
import OpenPNM.Utilities.transformations as tr


def voronoi(network, geometry, offset, **kwargs):
    r"""
    Offset the throat vertices effectively erroding the facet by the offset distance
    supplied.
    """
    Nt = geometry.num_throats()
    area = sp.ndarray(Nt)
    perimeter = sp.ndarray(Nt)
    offset_verts = sp.ndarray(Nt, dtype=object)
    offset_error = sp.ndarray(Nt)
    throat_COM = sp.ndarray([Nt, 3])
    for i in range(Nt):
        offset_rand = (sp.random.rand(1)-0.5)*offset
        throat_verts = geometry['throat.vertices'][i]
        throat_normal = geometry['throat.normal'][i]
        area[i], perimeter[i], offset_verts[i], throat_COM[i], offset_error[i] = \
            vo.get_throat_geom(throat_verts, throat_normal, offset)

    for i in range(Nt):
        if offset_error[i] > 0 and len(offset_verts[i]) > 0:
            offset_verts[i] = []
    # Properties that depend on the offset vertices are the area, perimeter and the
    # centroid or COM. To speed things up we could save them all now rather than
    # processing them individually.
    if kwargs['set_dependent'] is True:
        geometry['throat.area'] = area
        geometry['throat.perimeter'] = perimeter
        geometry['throat.centroid'] = throat_COM

    return offset_verts


def distance_transform(network, geometry, offset, **kwargs):
    r"""
    Use the Voronoi vertices and perform image analysis to obtain throat properties
    """

    import math
    import numpy as np
    from skimage.morphology import convex_hull_image
    from skimage.measure import regionprops
    from scipy import ndimage

    Nt = geometry.num_throats()
    area = sp.zeros(Nt)
    perimeter = sp.zeros(Nt)
    centroid = sp.zeros([Nt, 3])
    incentre = sp.zeros([Nt, 3])
    inradius = sp.zeros(Nt)
    equiv_diameter = sp.zeros(Nt)
    eroded_verts = sp.ndarray(Nt, dtype=object)

    res = 200
    vertices = geometry['throat.vertices']
    normals = geometry['throat.normal']
    z_axis = [0, 0, 1]

    for i in range(Nt):
        # For boundaries some facets will already be aligned with the axis - if this
        # is the case a rotation is unnecessary and could also cause problems
        angle = tr.angle_between_vectors(normals[i], z_axis)
        if angle == 0.0 or angle == np.pi:
            # We are already aligned
            rotate_facet = False
            facet = vertices[i]
        else:
            rotate_facet = True
            M = tr.rotation_matrix(tr.angle_between_vectors(normals[i], z_axis),
                                   tr.vector_product(normals[i], z_axis))
            facet = np.dot(vertices[i], M[:3, :3].T)
        x = facet[:, 0]
        y = facet[:, 1]
        z = facet[:, 2]
        # Get points in 2d for image analysis
        pts = np.column_stack((x, y))
        # Translate points so min sits at the origin
        translation = [pts[:, 0].min(), pts[:, 1].min()]
        pts -= translation
        order = np.int(math.ceil(-np.log10(np.max(pts))))
        # Normalise and scale the points so that largest span equals the resolution
        # to save on memory and create clear image"
        max_factor = np.max([pts[:, 0].max(), pts[:, 1].max()])
        f = res/max_factor
        # Scale the offset and define a circular structuring element with radius
        r = f*offset
        # Only proceed if r is less than half the span of the image"
        if r <= res/2:
            pts *= f
            minp1 = pts[:, 0].min()
            minp2 = pts[:, 1].min()
            maxp1 = pts[:, 0].max()
            maxp2 = pts[:, 1].max()
            img = np.zeros([np.int(math.ceil(maxp1-minp1)+1),
                            np.int(math.ceil(maxp2-minp2)+1)])
            int_pts = np.around(pts, 0).astype(int)
            for pt in int_pts:
                img[pt[0]][pt[1]] = 1
            # Pad with zeros all the way around the edges
            img_pad = np.zeros([np.shape(img)[0] + 2, np.shape(img)[1] + 2])
            img_pad[1:np.shape(img)[0]+1, 1:np.shape(img)[1]+1] = img

            # All points should lie on this plane but could be some rounding errors
            # so use the order parameter
            z_plane = sp.unique(np.around(z, order+2))
            if len(z_plane) > 1:
                print('Rotation for image analysis failed')
            "Fill in the convex hull polygon"
            convhullimg = convex_hull_image(img_pad)
            # Perform a Distance Transform and black out points less than r to create
            # binary erosion. This is faster than performing an erosion and dt can
            # also be used later to find incircle"
            eroded = ndimage.distance_transform_edt(convhullimg)
            eroded[eroded <= r] = 0
            eroded[eroded > r] = 1
            # If we are left with less than 3 non-zero points then the throat is
            # fully occluded
            if np.sum(eroded) >= 3:
                # Do some image analysis to extract the key properties
                regions = regionprops(eroded[1:np.shape(img)[0]+1,
                                             1:np.shape(img)[1]+1].astype(int))
                # Change this to cope with genuine multi-region throats
                if len(regions) == 1:
                    for props in regions:
                        x0, y0 = props.centroid
                        equiv_diameter[i] = props.equivalent_diameter
                        area[i] = props.area
                        perimeter[i] = props.perimeter
                        coords = props.coords
                    # Undo the translation, scaling and truncation on the centroid
                    centroid2d = [x0, y0]/f
                    centroid2d += (translation)
                    centroid3d = np.concatenate((centroid2d, z_plane))
                    # Distance transform the eroded facet to find the incentre and
                    # inradius
                    dt = ndimage.distance_transform_edt(eroded)
                    inx0, iny0 = \
                        np.asarray(np.unravel_index(dt.argmax(), dt.shape)) \
                          .astype(float)
                    incentre2d = [inx0, iny0]
                    # Undo the translation, scaling and truncation on the incentre
                    incentre2d /= f
                    incentre2d += (translation)
                    incentre3d = np.concatenate((incentre2d, z_plane))
                    # The offset vertices will be those in the coords that are
                    # closest to the originals"
                    offset_verts = []
                    for pt in int_pts:
                        vert = np.argmin(np.sum(np.square(coords-pt), axis=1))
                        if vert not in offset_verts:
                            offset_verts.append(vert)
                    # If we are left with less than 3 different vertices then the
                    # throat is fully occluded as we can't make a shape with
                    # non-zero area
                    if len(offset_verts) >= 3:
                        offset_coords = coords[offset_verts].astype(float)
                        # Undo the translation, scaling and truncation on the
                        # offset_verts
                        offset_coords /= f
                        offset_coords_3d = \
                            np.vstack((offset_coords[:, 0]+translation[0],
                                       offset_coords[:, 1]+translation[1],
                                       np.ones(len(offset_verts))*z_plane)).T

                        # Get matrix to un-rotate the co-ordinates back to the
                        # original orientation if we rotated in the first place
                        if rotate_facet:
                            MI = tr.inverse_matrix(M)
                            # Unrotate the offset coordinates
                            incentre[i] = np.dot(incentre3d, MI[:3, :3].T)
                            centroid[i] = np.dot(centroid3d, MI[:3, :3].T)
                            eroded_verts[i] = np.dot(offset_coords_3d, MI[:3, :3].T)

                        else:
                            incentre[i] = incentre3d
                            centroid[i] = centroid3d
                            eroded_verts[i] = offset_coords_3d

                        inradius[i] = dt.max()
                        # Undo scaling on other parameters
                        area[i] /= f*f
                        perimeter[i] /= f
                        equiv_diameter[i] /= f
                        inradius[i] /= f
                    else:
                        area[i] = 0
                        perimeter[i] = 0
                        equiv_diameter[i] = 0

    if kwargs['set_dependent'] is True:
        geometry['throat.area'] = area
        geometry['throat.perimeter'] = perimeter
        geometry['throat.centroid'] = centroid
        geometry['throat.diameter'] = equiv_diameter
        geometry['throat.indiameter'] = inradius*2
        geometry['throat.incentre'] = incentre

    return eroded_verts
