# -*- coding: utf-8 -*-
"""
===============================================================================
Delaunay: Generate random networks based on Delaunay Tessellations
===============================================================================

"""
import sys
import scipy as sp
import numpy as np
import OpenPNM.Utilities.vertexops as vo
import scipy.sparse as sprs
import scipy.spatial as sptl
import scipy.ndimage as spim
from scipy.spatial import Voronoi
from OpenPNM.Network import GenericNetwork
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class Delaunay(GenericNetwork):
    r"""
    This class contains the methods for creating a *Delaunay* network topology
    based connecting pores with a Delaunay tessellation.

    Parameters
    ----------
    name : string
        A unique name for the network
    domain_size : list of floats, [Lx,Ly,Lz]
        Bounding cube for internal pore positions
    num_pores : int
        Number of pores to place randomly within domain
    prob : 3D float array
        Values should be between 0 and 1 as determines probability of point
        with relative domain coordinates being kept. Array does not have to
        be same size as domain because positions are re-scaled
    base_points : [Np,3] float array
        coordinates to use instead of random generation
    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.Delaunay(num_pores=100,
    ...                               domain_size=[0.0001, 0.0001, 0.0001])
    >>> pn.num_pores()
    100

    """

    def __init__(self, num_pores=None, domain_size=None, prob=None,
                 base_points=None, **kwargs):
        """
        Create Delauny network object
        """
        super().__init__(**kwargs)
        self.generate(num_pores, domain_size, prob, base_points)

    def generate(self, num_pores, domain_size, prob, base_points):
        r"""
        Method to trigger the generation of the network
        """
        logger.info('Start of network topology generation')
        self._generate_setup(num_pores, domain_size, base_points)
        if base_points is not None:
            try:
                dim = sp.shape(base_points)[1]
                if dim != 3:
                    raise Exception('base points must be 3D')
            except:
                raise Exception('base points must be 3D')
            self['pore.coords'] = base_points
        else:
            self._generate_pores(prob)
        self._generate_throats()
        logger.debug('Network generation complete')

    def _generate_setup(self, num_pores, domain_size, base_points):
        r"""
        Perform applicable preliminary checks and calculations required for
        generation
        """
        logger.debug('generate_setup: Perform preliminary calculations and checks')
        if domain_size is None:
            raise Exception('domain_size must always be specified')
        if num_pores is None and base_points is None:
            raise Exception('num_pores or base_points must be specified')
        elif num_pores is None and base_points is not None:
            num_pores = len(base_points)
        elif num_pores is not None and base_points is not None:
            logger.warning('both num_pores and base_points arguments given' +
                           ' num_pores over-written')
            num_pores = len(base_points)

        self._Lx = domain_size[0]
        self._Ly = domain_size[1]
        self._Lz = domain_size[2]
        self._Np = num_pores
        r"""
        TODO: Fix this, btype should be received as an argument
        """
        self._btype = [0, 0, 0]

    def _generate_pores(self, prob=None):
        r"""
        Generate the pores with numbering scheme.
        """
        logger.info('Place randomly located pores in the domain')
        if prob is not None:
            coords = []
            i = 0
            while i < self._Np:
                coord = np.random.rand(3)
                [indx, indy, indz] = np.floor(coord*np.shape(prob)).astype(int)
                p = prob[indx][indy][indz]
                if np.random.rand(1) <= p:
                    coords.append(coord)
                    i += 1
            coords = np.asarray(coords)
        else:
            coords = np.random.random([self._Np, 3])

        coords *= np.array([self._Lx, self._Ly, self._Lz])
        self['pore.coords'] = coords

    def _generate_throats(self):
        r"""
        Generate the throats connections
        """
        logger.info('Define connections between pores')
        pts = self['pore.coords']
        Np = len(pts)
        # Generate 6 dummy domains to pad onto each face of real domain This
        # prevents surface pores from making long range connections to each other

        x, y, z = self['pore.coords'].T
        if x.max() > self._Lx:
            Lx = x.max()*1.05
        else:
            Lx = self._Lx
        if y.max() > self._Ly:
            Ly = y.max()*1.05
        else:
            Ly = self._Ly
        if z.max() > self._Lz:
            Lz = z.max()*1.05
        else:
            Lz = self._Lz

        # Reflect in X = Lx and 0
        Pxp = pts.copy()
        Pxp[:, 0] = 2*Lx-Pxp[:, 0]
        Pxm = pts.copy()
        Pxm[:, 0] = Pxm[:, 0]*(-1)
        # Reflect in Y = Ly and 0
        Pyp = pts.copy()
        Pyp[:, 1] = 2*Ly-Pxp[:, 1]
        Pym = pts.copy()
        Pym[:, 1] = Pxm[:, 1]*(-1)
        # Reflect in Z = Lz and 0
        Pzp = pts.copy()
        Pzp[:, 2] = 2*Lz-Pxp[:, 2]
        Pzm = pts.copy()
        Pzm[:, 2] = Pxm[:, 2]*(-1)
        # Add dummy domains to real domain
        # Order important for boundary logic
        pts = np.vstack((pts, Pxp, Pxm, Pyp, Pym, Pzp, Pzm))
        # Perform tessellation
        logger.debug('Beginning tessellation')
        Tri = sptl.Delaunay(pts)
        logger.debug('Converting tessellation to adjacency matrix')
        adjmat = sprs.lil_matrix((Np, Np), dtype=int)
        for i in sp.arange(0, sp.shape(Tri.simplices)[0]):
            # Keep only simplices that are fully in real domain
            # this used to be vectorize, but it stopped working...change in scipy?
            for j in Tri.simplices[i]:
                if j < Np:
                    adjmat[j, Tri.simplices[i][Tri.simplices[i] < Np]] = 1
        # Remove duplicate (lower triangle) and self connections (diagonal)
        # and convert to coo
        adjmat = sprs.triu(adjmat, k=1, format="coo")
        logger.debug('Conversion to adjacency matrix complete')
        self['throat.conns'] = sp.vstack((adjmat.row, adjmat.col)).T
        self['pore.all'] = np.ones(len(self['pore.coords']), dtype=bool)
        self['throat.all'] = np.ones(len(self['throat.conns']), dtype=bool)

        # Do Voronoi diagram - creating voronoi polyhedra around each pore and save
        # vertex information
        self._vor = Voronoi(pts)
        all_vert_index = sp.ndarray(Np, dtype=object)
        for i, polygon in enumerate(self._vor.point_region[0:Np]):
            if -1 not in self._vor.regions[polygon]:
                all_vert_index[i] = \
                    dict(zip(self._vor.regions[polygon],
                             self._vor.vertices[self._vor.regions[polygon]]))

        # Add throat vertices by looking up vor.ridge_dict
        throat_verts = sp.ndarray(len(self['throat.conns']), dtype=object)
        for i, (p1, p2) in enumerate(self['throat.conns']):
            try:
                throat_verts[i] = \
                    dict(zip(self._vor.ridge_dict[(p1, p2)],
                             self._vor.vertices[self._vor.ridge_dict[(p1, p2)]]))
            except KeyError:
                try:
                    throat_verts[i] = \
                        dict(zip(self._vor.ridge_dict[(p2, p1)],
                                 self._vor.vertices[self._vor.ridge_dict[(p2, p1)]]))
                except KeyError:
                    logger.error('Throat Pair Not Found in Voronoi Ridge Dictionary')

        self['pore.vert_index'] = all_vert_index
        self['throat.vert_index'] = throat_verts
        logger.debug(sys._getframe().f_code.co_name + ': End of method')

    def add_boundaries(self):

        r"""
        This method identifies pores in the original Voronoi object that straddle a
        boundary imposed by the reflection. The pore inside the original set of pores
        (with index 0 - Np) is identified and the coordinates are saved. The vertices
        making up the boundary throat are retrieved from the ridge_dict values and
        these are used to identify which boundary the throat sits at.
        A new pore and new connection is created with coordinates lying on the
        boundary plane.
        N.B This method will only work properly if the original network remains
            unaltered i.e. not trimmed or extended
            This preserves the connection between pore index on the network object
            and the Voronoi object
            The point of using this method is so that the throat vertices created by
            the Voronoi object are preserved

        This method will create boundary pores at the centre of the voronoi faces
        that align with the outer planes of the domain.
        The original pores in the domain are labelled internal and the boundary pores
        are labelled external

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.Delaunay(num_pores=100,
        ...                               domain_size=[0.0001,0.0001,0.0001])
        >>> pn.add_boundaries()
        >>> pn.num_pores('boundary') > 0
        True
        """

        bound_conns = []
        bound_coords = []
        bound_vert_index = []
        throat_vert_index = []
        # Find boundary extent
        [x_min, x_max, y_min, y_max, z_min, z_max] = \
            vo.vertex_dimension(self, self.pores(), parm='minmax')
        min_point = np.around(np.array([x_min, y_min, z_min]), 10)
        max_point = np.around(np.array([x_max, y_max, z_max]), 10)
        Np = self.num_pores()
        Nt = self.num_throats()
        new_throat_count = 0
        # ridge_dict contains a dictionary where the key is a set of 2 neighbouring
        # pores and the value is the vertex indices that form the throat or ridge
        # between them
        for p, v in self._vor.ridge_dict.items():
            # If the vertex with index -1 is contained in list then the ridge is
            # unbounded - ignore these
            if np.all(np.asarray(v) >= 0):
                # Boundary throats will be those connecting one pore inside the
                # original set and one out
                if (p[0] in range(Np) and p[1] not in range(Np)) or \
                        (p[0] not in range(Np) and p[1] in range(Np)):
                    # The dictionary key is not in numerical order so find the pore
                    # index inside
                    if p[0] in range(Np):
                        my_pore = p[0]
                    else:
                        my_pore = p[1]
                    my_pore_coord = self["pore.coords"][my_pore]
                    new_pore_coord = my_pore_coord.copy()
                    # Rounding necessary here to identify the plane as Voronoi can
                    # have 1e-17 and smaller errors
                    throat_verts = np.around(self._vor.vertices[v], 10)
                    # Find which plane we are aligned with (if any) and align
                    # new_pore with throat plane
                    if len(np.unique(throat_verts[:, 0])) == 1:
                        new_pore_coord[0] = np.unique(throat_verts[:, 0])
                    elif len(np.unique(throat_verts[:, 1])) == 1:
                        new_pore_coord[1] = np.unique(throat_verts[:, 1])
                    elif len(np.unique(throat_verts[:, 2])) == 1:
                        new_pore_coord[2] = np.unique(throat_verts[:, 2])
                    else:
                        new_pore_coord = np.mean(throat_verts, axis=0)
                        pass
                    bound_coords.append(new_pore_coord)
                    bound_conns.append(np.array([my_pore, new_throat_count + Np]))
                    bound_vert_index.append(dict(zip(v, throat_verts)))
                    throat_vert_index.append(dict(zip(v, throat_verts)))
                    new_throat_count += 1

        # Add new pores and connections
        self.extend(pore_coords=bound_coords, throat_conns=bound_conns)
        # Record new number of pores
        Mp = self.num_pores()
        Mt = self.num_throats()
        new_pore_ids = np.arange(Np, Mp)
        new_throat_ids = np.arange(Nt, Mt)
        # Identify which boundary the pore sits on
        front = self.pores()[self['pore.coords'][:, 0] == min_point[0]]
        back = self.pores()[self['pore.coords'][:, 0] == max_point[0]]
        left = self.pores()[self['pore.coords'][:, 1] == min_point[1]]
        right = self.pores()[self['pore.coords'][:, 1] == max_point[1]]
        bottom = self.pores()[self['pore.coords'][:, 2] == min_point[2]]
        top = self.pores()[self['pore.coords'][:, 2] == max_point[2]]
        if len(top) == 0:
            top = self.pores()[self['pore.coords'][:, 2] ==
                               np.asarray(bound_coords)[:, 2].max()]
        # Assign labels
        self['pore.boundary'] = False
        self['pore.boundary'][new_pore_ids] = True
        self['throat.boundary'] = False
        self['throat.boundary'][new_throat_ids] = True
        self['pore.right_boundary'] = False
        self['pore.left_boundary'] = False
        self['pore.front_boundary'] = False
        self['pore.back_boundary'] = False
        self['pore.top_boundary'] = False
        self['pore.bottom_boundary'] = False
        self['pore.right_boundary'][right] = True
        self['pore.left_boundary'][left] = True
        self['pore.front_boundary'][front] = True
        self['pore.back_boundary'][back] = True
        self['pore.top_boundary'][top] = True
        self['pore.bottom_boundary'][bottom] = True
        # Save the throat verts
        self["pore.vert_index"][new_pore_ids] = bound_vert_index
        self["throat.vert_index"][new_throat_ids] = throat_vert_index

    def domain_length(self, face_1, face_2):
        r"""
        Returns the distance between two faces
        No coplanar checking this is done in vertex_dimension
        """
        L = vo.vertex_dimension(self, face_1, face_2, parm='length')
        return L

    def domain_area(self, face):
        r"""
        Returns the area of a face
        No coplanar checking this is done in vertex_dimension
        """
        A = vo.vertex_dimension(self, face, parm='area')

        return A

    def _export_vor_fibres(self):
        r"""
        Run through the throat vertices, compute the convex hull order and save
        the vertices and ordered faces in a pickle dictionary to be used in
        blender
        """
        import pickle as pickle
        Indices = []
        for t in self.throats():
            indices = list(self["throat.vert_index"][t].keys())
            verts = self._vor.vertices[indices]
            # Need to order the indices in convex hull order
            # Compute the standard deviation in all coordinates and eliminate
            # the axis with the smallest to make 2d
            stds = [np.std(verts[:, 0]), np.std(verts[:, 1]), np.std(verts[:, 2])]
            if np.argmin(stds) == 0:
                verts2d = np.vstack((verts[:, 1], verts[:, 2])).T
            elif np.argmin(stds) == 1:
                verts2d = np.vstack((verts[:, 0], verts[:, 2])).T
            else:
                verts2d = np.vstack((verts[:, 0], verts[:, 1])).T
            # 2d convexhull returns vertices in hull order
            hull2d = sptl.ConvexHull(verts2d, qhull_options='QJ Pp')
            # Re-order the vertices and save as list (blender likes them as lists)
            Indices.append(np.asarray(indices)[hull2d.vertices].tolist())
        # Create dictionary to pickle
        data = {}
        data["Verts"] = self._vor.vertices
        data["Indices"] = Indices
        pickle.dump(data, open("fibres.p", "wb"))
