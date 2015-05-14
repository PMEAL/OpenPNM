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
from scipy import stats as st
from scipy.special import cbrt
from OpenPNM.Network import GenericNetwork
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class Delaunay(GenericNetwork):
    r"""
    This class contains the methods for creating a *Delaunay* network topology
    based connecting pores with a Delaunay tessellation.

    To invoke the actual generation it is necessary to run the `generate` method.

    Parameters
    ----------
    name : string
        A unique name for the network

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.Delaunay(num_pores=100,
    ...                               domain_size=[0.0001, 0.0001, 0.0001])
    >>> pn.num_pores()
    100

    """

    def __init__(self, num_pores=None, domain_size=None, **kwargs):
        """
        Create Delauny network object
        """
        super().__init__(**kwargs)
        if (num_pores and domain_size) is None:
            num_pores = 1
            domain_size = [1.0, 1.0, 1.0]
        else:
            self.generate(num_pores, domain_size)

    def generate(self, num_pores, domain_size):
        r"""
        Method to trigger the generation of the network

        Parameters
        ----------
        domain_size : list of floats, [Lx,Ly,Lz]
            Bounding cube for internal pore positions
        num_pores : int
            Number of pores to place randomly within domain

        """
        logger.info('Start of network topology generation')
        self._generate_setup(num_pores, domain_size)
        self._generate_pores()
        self._generate_throats()
        logger.debug('Network generation complete')

    def _generate_setup(self, num_pores, domain_size):
        r"""
        Perform applicable preliminary checks and calculations required for
        generation
        """
        logger.debug('generate_setup: Perform preliminary calculations')
        if domain_size is not None and num_pores is not None:
            self._Lx = domain_size[0]
            self._Ly = domain_size[1]
            self._Lz = domain_size[2]
            self._Np = num_pores
            r"""
            TODO: Fix this, btype should be received as an argument
            """
            self._btype = [0, 0, 0]
        else:
            logger.error('domain_size and num_pores must be specified')
            raise Exception('domain_size and num_pores must be specified')

    def _generate_pores(self):
        r"""
        Generate the pores with numbering scheme.
        """
        logger.info('Place randomly located pores in the domain')
        coords = np.zeros([self._Np, 3])
        # Reject points close to boundaries - if False there will be slightly more
        rejection = [False, False, True]
        for j in range(3):
            i = 0
            while i < self._Np:
                coord = np.random.uniform(0, 1, 1)
                if self._reject(coord) == rejection[j]:
                    coords[i][j] = coord
                    i += 1
        coords *= np.array([self._Lx, self._Ly, self._Lz])

        self['pore.coords'] = coords

    def _prob_func(self, m):
        a = 35
        b = 0.2
        p = ((m**a) + ((1-m)**a) + (2*b))/(1 + (2*b))
        return p

    def _reject(self, point):

        P = self._prob_func(point)
        nrand = np.random.uniform(0, 1, 1)
        # Place more points at the sides of the domain and fewer at the
        # top and bottom
        if P < nrand:
            rejection = True
        else:
            rejection = False

        return rejection

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
                    print('Throat Pair Not Found in Voronoi Ridge Dictionary')

        self['pore.vert_index'] = all_vert_index
        self['throat.vert_index'] = throat_verts
        logger.debug(sys._getframe().f_code.co_name + ': End of method')

    def _add_labels(self):
        r"""
        Deprecated if using add_boundaries()
        This finds surface pores simply by proximity to the domain boundaries.
        A better approach is necessary
        """
        coords = self['pore.coords']
        self['pore.front'] = coords[:, 0] < (0.1*self._Lx)
        self['pore.back'] = coords[:, 0] > (0.9*self._Lx)
        self['pore.left'] = coords[:, 1] < (0.1*self._Ly)
        self['pore.right'] = coords[:, 1] > (0.9*self._Ly)
        self['pore.bottom'] = coords[:, 2] < (0.1*self._Lz)
        self['pore.top'] = coords[:, 2] > (0.9*self._Lz)
        bnds = self.pores(labels=['front', 'back', 'left', 'right', 'bottom', 'top'])
        self['pore.boundary'] = False
        self['pore.boundary'] = bnds

    def _add_boundaries(self):
        r"""
        This is an alternative means of adding boundaries
        """
        logger.info('add_boundaries: start of method')

        import scipy.spatial as sptl
        import scipy.sparse as sprs
        Lx = self._Lx
        Ly = self._Ly
        Lz = self._Lz
        Np = self.num_pores()
        btype = self._btype
        boffset = 0.05

        # Translate internal pores to each face of domain
        poffset = np.zeros((7, 3))
        poffset[[2, 5], 0] = [-Lx, Lx]
        poffset[[3, 4], 1] = [-Ly, Ly]
        poffset[[1, 6], 2] = [-Lz, Lz]
        pcoords = pcoords0 = self['pore.coords']
        for i in np.r_[1:7]:
            pcoords = np.concatenate((pcoords, pcoords0 + poffset[i, :]), axis=0)

        # Use some twisted logic to get bval list of + for boundary and -
        # for periodic faces
        bval = [0, 1, 2, 3, 4, 5, 6] * \
            (np.array([0, btype[2], btype[0], btype[1],
                       btype[1], btype[0], btype[2]])*-2+1)
        ptype = np.zeros((Np,), dtype=int)
        for i in np.r_[1:7]:
            ptype = \
                np.concatenate((ptype, np.ones((Np,), dtype=int)*bval[i]), axis=0)

        # pnum contains the internal ID number of the boundary pores
        # for connecting periodic points
        pnum = self.pores()
        pnum = np.tile(pnum, 7)

        Tri = sptl.Delaunay(pcoords)
        adjmat = \
            sprs.lil_matrix((np.shape(pcoords)[0], np.shape(pcoords)[0]), dtype=int)
        for i in np.arange(0, np.shape(Tri.simplices)[0]):
            # Keep only simplices that are fully in real domain
            adjmat[Tri.simplices[i], Tri.simplices[i]] = 1
        adjmat = sprs.triu(adjmat, k=1, format='lil')
        for i in np.arange(0, Np):
            # Add periodic throats to the netowrk (if any)
            tpore2 = pnum[adjmat.rows[i]][ptype[adjmat.rows[i]] < 0]
            tpore1 = np.ones_like(tpore2, dtype=int) * i
            conns = self['throat.conns']
            conns = np.concatenate((conns, np.vstack((tpore1, tpore2)).T), axis=0)
            # Add boundary pores and throats to the network
            newporetyps = np.unique(ptype[adjmat.rows[i]][ptype[adjmat.rows[i]] > 0])
            newporenums = \
                np.r_[self.num_pores():self.num_pores()+np.size(newporetyps)]
            tpore2 = newporenums
            tpore1 = np.ones_like(tpore2, dtype=int) * i
            conns = np.concatenate((conns, np.vstack((tpore1, tpore2)).T), axis=0)
            self['throat.conns'] = conns
            bcoords = np.zeros((7, 3), dtype=float)
            coords = self['pore.coords']
            bcoords[1, :] = [coords[i, 0], coords[i, 1], 0-Lz*boffset]
            bcoords[2, :] = [0-Lx*boffset, coords[i, 1], coords[i, 2]]
            bcoords[3, :] = [coords[i, 0], -Ly*boffset, coords[i, 2]]
            bcoords[4, :] = [coords[i, 0], Ly+Ly*boffset, coords[i, 2]]
            bcoords[5, :] = [Lx+Lx*boffset, coords[i, 1], coords[i, 2]]
            bcoords[6, :] = [coords[i, 0], coords[i, 1], Lz+Lz*boffset]
            newporecoords = bcoords[newporetyps, :]
            coords = np.concatenate((coords, newporecoords), axis=0)
            self['pore.coords'] = coords
        # Reset number of pores and throats (easier than tracking it)
        nums = np.r_[0:np.shape(coords)[0]]
        self['pore.numbering'] = nums
        self['pore.numbering'] = np.ones((nums[-1]+1,), dtype=bool)
        nums = np.r_[0:np.shape(conns)[0]]
        self['throat.numbering'] = nums
        self['throat.numbering'] = np.ones((nums[-1]+1,), dtype=bool)
        logger.debug('add_boundaries: end of method')

    def _add_boundaries_old(self):
        logger.info('add_boundaries_old: Start of method')

        self.add_opposing_boundaries(btype=[2, 5])
        self.add_opposing_boundaries(btype=[3, 4])
        self.add_opposing_boundaries(btype=[1, 6])

    def _add_opposing_boundaries(self, btype=[1, 6]):
        r"""
        btype indicates which two boundaries are being added by type
        """
        logger.info('add_opposing_boundaries: start of method')

        if btype == [2, 5]:
            D = 0
            W = 1
            H = 2
        elif btype == [3, 4]:
            D = 1
            W = 0
            H = 2
        elif btype == [1, 6]:
            D = 2
            W = 1
            H = 0

        Lx = self.domain_size[D]
        Ly = self.domain_size[W]
        Lz = self.domain_size[H]
        # Rotate pore coordinates (use only internal pores)
        pnum = self._net.pore_data['numbering'][self._net.pore_data['type'] == 0]
        pcoords = np.zeros_like(self._net.pore_data['coords'][pnum, :])
        pcoords[:, 0] = self._net.pore_data['coords'][pnum, D]
        pcoords[:, 1] = self._net.pore_data['coords'][pnum, W]
        pcoords[:, 2] = self._net.pore_data['coords'][pnum, H]

        # Determine dimensions of image from dimensions of domain
        f = 100  # minimum image dimension
        im_dim = [0, 0, 0]
        im_dim[0] = np.floor(f*Lx/np.min([Lx, Ly, Lz]))
        im_dim[1] = np.floor(f*Ly/np.min([Lx, Ly, Lz]))
        im_dim[2] = np.floor(f*Lz/np.min([Lx, Ly, Lz]))
        im_dim = np.array(im_dim, dtype=int)

        # Convert pore coordinates into image subscripts
        im_subs = np.zeros_like(pcoords, dtype=int)
        im_subs[:, 0] = pcoords[:, 0]*im_dim[0]/Lx
        im_subs[:, 1] = pcoords[:, 1]*im_dim[1]/Ly
        im_subs[:, 2] = pcoords[:, 2]*im_dim[2]/Lz
        # Find linear indices of each pore in the new image
        im_inds = np.ravel_multi_index((im_subs[:, 0], im_subs[:, 1], im_subs[:, 2]),
                                       dims=(im_dim), order='F')

        # Generate 3D image of points (place pore numbers at each site for use later)
        img = np.zeros(im_dim, dtype=int)
        img[im_subs[:, 0], im_subs[:, 1], im_subs[:, 2]] = pnum

        # Perform distance transform on points and also get 'indicies' of each point
        img_dt, ind_dt = spim.distance_transform_edt(img == 0)

        # Project all* internal points to x face
        # *Note that it's possible/likely that mutliple internal points map to the
        # same boundary point
        img_bd0 = np.zeros([im_dim[1], im_dim[2]], dtype=int)
        img_bd1 = np.zeros([im_dim[1], im_dim[2]], dtype=int)
        img_bd0[im_subs[:, 1], im_subs[:, 2]] = im_inds
        img_bd1[im_subs[:, 1], im_subs[:, 2]] = im_inds

        # Create 2D array of distance transform indices for 0 and end faces
        dt_D0 = ind_dt[0, 0, :, :]*(img_bd0 > 0)  # 0 face
        dt_D1 = ind_dt[0, -1, :, :]*(img_bd1 > 0)  # end face

        # Create a 2D mask containing x coordinates of internal points
        # and -1 elsewhere
        img_D0 = -np.ones([im_dim[1], im_dim[2]], dtype=int)
        img_D1 = -np.ones([im_dim[1], im_dim[2]], dtype=int)
        img_D0[im_subs[:, 1], im_subs[:, 2]] = im_subs[:, 0]
        img_D1[im_subs[:, 1], im_subs[:, 2]] = im_subs[:, 0]

        # Find where x value of internal points corresponds to x value of
        # distance transform indices
        img_bd0 = (img_D0 == dt_D0)*img_bd0
        img_bd1 = (img_D1 == dt_D1)*img_bd1

        # Convert boundary sites to linear indices
        inds_bd0 = img_bd0[np.nonzero(img_bd0)]
        inds_bd1 = img_bd1[np.nonzero(img_bd1)]

        # Use linear indices to find pore ID nums
        nums_bd0 = img[np.unravel_index(inds_bd0, dims=(im_dim), order='F')]
        nums_bd1 = img[np.unravel_index(inds_bd1, dims=(im_dim), order='F')]
        nums_bd = np.append(nums_bd0, nums_bd1)
        types_bd = np.append(np.zeros_like(nums_bd0), np.ones_like(nums_bd1))

        # Add new boundary pores and throats to the network
        # Get all pores including previously added boundaries
        Np = self._net.num_pores()
        bp_numbering = np.r_[Np:Np+np.size(nums_bd)]
        bp_type = (types_bd == 0)*btype[0] + (types_bd == 1)*btype[1]
        bp_coords = np.zeros([np.size(nums_bd), 3])
        bp_coords[types_bd == 0, D] = np.zeros_like(nums_bd0) - .0001
        bp_coords[types_bd == 0, W] = pcoords[nums_bd0, 1]
        bp_coords[types_bd == 0, H] = pcoords[nums_bd0, 2]
        bp_coords[types_bd == 1, D] = np.ones_like(nums_bd1)*Lx+0.0001
        bp_coords[types_bd == 1, W] = pcoords[nums_bd1, 1]
        bp_coords[types_bd == 1, H] = pcoords[nums_bd1, 2]
        self._net.pore_data['numbering'] = \
            np.append(self._net.pore_data['numbering'], bp_numbering)
        self._net.pore_data['type'] = np.append(self._net.pore_data['type'], bp_type)
        self._net.pore_data['coords'] = \
            np.concatenate((self._net.pore_data['coords'], bp_coords))
        Nt = self._net.num_throats()
        bt_numbering = np.r_[Nt:Nt + np.size(nums_bd)]
        bt_type = np.ones(np.size(nums_bd), dtype=int)*2
        bt_connections = np.zeros([np.size(nums_bd), 2], dtype=int)
        bt_connections[:, 0] = nums_bd
        bt_connections[:, 1] = bp_numbering
        self._net.throat_data['numbering'] = \
            np.append(self._net.throat_data['numbering'], bt_numbering)
        self._net.throat_data['type'] = \
            np.append(self._net.throat_data['type'], bt_type)
        self._net.throat_data['conns'] = \
            np.concatenate((self._net.throat_data['conns'], bt_connections))

    def domain_size(self, dimension=''):
        r"""
        This is a simple way to find the domain sizes.
        N.B
        Will not work with saved and loaded networks
        """
        if dimension == 'front' or dimension == 'back':
            return self._Ly*self._Lz
        if dimension == 'left' or dimension == 'right':
            return self._Lx*self._Lz
        if dimension == 'top' or dimension == 'bottom':
            return self._Lx*self._Ly
        if dimension == 'volume':
            return self._Lx*self._Ly*self._Lz
        if dimension == 'height':
            return self._Lz
        if dimension == 'width':
            return self._Lx
        if dimension == 'depth':
            return self._Ly

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
                        new_pore_coord = throat_verts.mean()
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
        # Assign labels
        self['pore.boundary'] = False
        self['pore.boundary'][new_pore_ids] = True
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

        Example
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.Delaunay(num_pores=100, domain_size=[3,2,1])
        >>> pn.add_boundaries()
        >>> B1 = pn.pores('left_boundary')
        >>> B2 = pn.pores('right_boundary')
        >>> pn.domain_length(B1,B2)
        2.0
        """
        L = vo.vertex_dimension(self, face_1, face_2, parm='length')
        return L

    def domain_area(self, face):
        r"""
        Returns the area of a face
        No coplanar checking this is done in vertex_dimension
        Example
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.Delaunay(num_pores=100, domain_size=[3,2,1])
        >>> pn.add_boundaries()
        >>> B1 = pn.pores('left_boundary')
        >>> B2 = pn.pores('right_boundary')
        >>> pn.domain_area(B1)
        3.0
        """
        A = vo.vertex_dimension(self, face, parm='area')

        return A

    def trim_occluded_throats(self):
        r"""
        After the offsetting routine throats with zero area have been fully occluded.
        Remove these from the network and also remove pores that are isolated
        """
        occluded_ts = list(self.throats()[self['throat.area'] == 0])
        if len(occluded_ts) > 0:
            self.trim(throats=occluded_ts)
        # Also get rid of isolated pores
        isolated_ps = self.check_network_health()['isolated_pores']
        if len(isolated_ps) > 0:
            self.trim(isolated_ps)
