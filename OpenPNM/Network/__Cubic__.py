# -*- coding: utf-8 -*-
"""
===============================================================================
Cubic: Generate lattice-like networks
===============================================================================

"""
import numpy as np
import scipy as sp
import scipy.spatial as sptl
import OpenPNM.Utilities.misc as misc
from OpenPNM.Network import tools
from OpenPNM.Network import GenericNetwork
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class Cubic(GenericNetwork):
    r"""
    This class generates a cubic network of the specified size and shape.
    Alternatively, an arbitrary domain shape defined by a supplied template.

    Parameters
    ----------
    name : string
        A unique name for the network

    shape : tuple of ints
        The (i,j,k) size and shape of the network.

    connectivity : int
        The number of connections to neighboring pores.  Connections are made
        symmetrically to any combination of face, edge or corners neighbors.

        Options are:

        - 6: Faces only
        - 8: Corners only
        - 12: Edges Only
        - 14: Faces and Corners
        - 18: Faces and Edges
        - 20: Edges and Corners
        - 26: Faces, Edges and Corners

    template : array of booleans
        An (i,j,k) array with True where the Network should be defined and
        False elsewhere. This approach is useful for creating networks of non-
        cuboid shape like spheres or cylinders, but still with a cubic lattice
        topology.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.Cubic(shape=[3,4,5])
    >>> pn.Np
    60

    It is also possible to create Networks with cubic connectivity but
    non-Cubic shape by provding an array with True values where the network
    should exist to the ``template`` argument. The example below produces a sphere:

    >>> img = sp.ones([11, 11, 11])
    >>> img[5, 5, 5] = 0
    >>> from scipy.ndimage import distance_transform_bf as dt
    >>> img = dt(img) < 5  # Create a sphere of True
    >>> pn = OpenPNM.Network.Cubic(template=img)
    >>> pn.Np
    485

    If random distributions of coordination number is desired, one option is
    to create a Cubic network with many connections and the trim some:

    >>> pn = OpenPNM.Network.Cubic(shape=[5, 5, 5], connectivity=26)
    >>> Nt_original = pn.Nt
    >>> mod = OpenPNM.Network.models.pore_topology.reduce_coordination
    >>> pn.add_model(propname='throat.to_drop', model=mod, z=10, mode='random')
    >>> pn.trim(throats=pn['throat.to_drop'])
    >>> pn.Nt < Nt_original
    True
    """
    def __init__(self, shape=None, template=None, spacing=[1, 1, 1],
                 connectivity=6, **kwargs):
        super().__init__(**kwargs)

        if shape is not None:
            arr = np.atleast_3d(np.empty(shape))
        elif template is not None:
            arr = sp.array(template, ndmin=3, dtype=bool)
        else:
            arr = np.atleast_3d(np.empty([1, 1, 1]))

        # Store original network shape
        self._shape = sp.shape(arr)
        # Store network spacing
        self._spacing = sp.ones(3)*sp.array(spacing, ndmin=1)

        points = np.array([i for i, v in np.ndenumerate(arr)], dtype=float)
        points += 0.5
        points *= spacing

        I = np.arange(arr.size).reshape(arr.shape)

        face_joints = [
            (I[:, :, :-1], I[:, :, 1:]),
            (I[:, :-1], I[:, 1:]),
            (I[:-1], I[1:]),
        ]

        corner_joints = [
            (I[:-1, :-1, :-1], I[1:, 1:, 1:]),
            (I[:-1, :-1, 1:], I[1:, 1:, :-1]),
            (I[:-1, 1:, :-1], I[1:, :-1, 1:]),
            (I[1:, :-1, :-1], I[:-1, 1:, 1:]),
        ]

        edge_joints = [
            (I[:, :-1, :-1], I[:, 1:, 1:]),
            (I[:, :-1, 1:], I[:, 1:, :-1]),
            (I[:-1, :, :-1], I[1:, :, 1:]),
            (I[1:, :, :-1], I[:-1, :, 1:]),
            (I[1:, 1:, :], I[:-1, :-1, :]),
            (I[1:, :-1, :], I[:-1, 1:, :]),
        ]

        if connectivity == 6:
            joints = face_joints
        elif connectivity == 8:
            joints = corner_joints
        elif connectivity == 12:
            joints = edge_joints
        elif connectivity == 14:
            joints = face_joints + corner_joints
        elif connectivity == 18:
            joints = face_joints + edge_joints
        elif connectivity == 20:
            joints = edge_joints + corner_joints
        elif connectivity == 26:
            joints = face_joints + corner_joints + edge_joints
        else:
            raise Exception('Invalid connectivity receieved. Must be 6, 8, 12, 14, '
                            '18, 20 or 26')

        I = np.arange(arr.size).reshape(arr.shape)
        tails, heads = [], []
        for T, H in joints:
            tails.extend(T.flat)
            heads.extend(H.flat)
        pairs = np.vstack([tails, heads]).T

        self['pore.coords'] = points
        self['throat.conns'] = pairs
        self['pore.all'] = np.ones(len(self['pore.coords']), dtype=bool)
        self['throat.all'] = np.ones(len(self['throat.conns']), dtype=bool)
        self['pore.index'] = sp.arange(0, len(self['pore.coords']))

        self._label_surfaces()

        # If an image was sent as 'template', then trim network to image shape
        if template is not None:
            self.trim(~arr.flatten())

    def _label_surfaces(self):
        r'''
        It applies the default surface labels for a cubic network
        '''
        x, y, z = self['pore.coords'].T
        labels = ['internal', 'front', 'back', 'left', 'right', 'bottom', 'top']
        for label in labels:
            if 'pore.'+label not in self.keys():
                self['pore.'+label] = False
        if 'pore.boundary' in self.keys():
            internal = -self['pore.boundary']
        else:
            internal = self['pore.all']
        self['pore.internal'] = internal
        self['pore.front'][x <= x.min()] = True
        self['pore.back'][x >= x.max()] = True
        self['pore.left'][y <= y.min()] = True
        self['pore.right'][y >= y.max()] = True
        self['pore.bottom'][z <= z.min()] = True
        self['pore.top'][z >= z.max()] = True

    def add_boundaries(self):
        r"""
        This method uses ``clone_pores`` to clone the surface pores (labeled
        'left','right', etc), then shifts them to the periphery of the domain,
        and gives them the label 'right_face', 'left_face', etc.
        """
        x, y, z = self['pore.coords'].T
        Lcx, Lcy, Lcz = self._spacing

        offset = {}
        offset['front'] = offset['left'] = offset['bottom'] = [0, 0, 0]
        offset['back'] = [Lcx*self._shape[0], 0, 0]
        offset['right'] = [0, Lcy*self._shape[1], 0]
        offset['top'] = [0, 0, Lcz*self._shape[2]]

        scale = {}
        scale['front'] = scale['back'] = [0, 1, 1]
        scale['left'] = scale['right'] = [1, 0, 1]
        scale['bottom'] = scale['top'] = [1, 1, 0]

        for label in ['front', 'back', 'left', 'right', 'bottom', 'top']:
            ps = self.pores(label)
            self.clone_pores(pores=ps, apply_label=[label+'_boundary', 'boundary'])
            # Translate cloned pores
            ind = self.pores(label+'_boundary')
            coords = self['pore.coords'][ind]
            coords = coords*scale[label] + offset[label]
            self['pore.coords'][ind] = coords

    def add_boundary_pores(self, pores, offset, apply_label='boundary'):
        r"""
        This method uses ``clone_pores`` to clone the input pores, then shifts
        them the specified amount and direction, then applies the given label.

        Parameters
        ----------
        pores : array_like
            List of pores to offset.  If no pores are specified, then it
            assumes that all surface pores are to be cloned.

        offset : 3 x 1 array
            The distance in vector form which the cloned boundary pores should
            be offset.  If no spacing is provided, then the spacing is inferred
            from the Network.

        apply_label : string
            This label is applied to the boundary pores.  Default is
            'boundary'.

        Examples
        --------
        >>> import OpenPNM as op
        >>> pn = op.Network.Cubic(shape=[5, 5, 5])
        >>> print(pn.Np)  # Confirm initial Network size
        125
        >>> Ps = pn.pores('top')  # Select pores on top face
        >>> pn.add_boundary_pores(pores=Ps, offset=[0, 0, 1])
        >>> print(pn.Np)  # Confirm addition of 25 new pores
        150
        >>> 'pore.boundary' in pn.labels()  # Default label is created
        True
        """
        # Parse the input pores
        Ps = sp.array(pores, ndmin=1)
        if Ps.dtype is bool:
            Ps = self.toindices(Ps)
        if sp.size(pores) == 0:  # Handle an empty array if given
            return sp.array([], dtype=sp.int64)
        # Clone the specifed pores
        self.clone_pores(pores=Ps)
        newPs = self.pores('pore.clone')
        del self['pore.clone']
        # Offset the cloned pores
        self['pore.coords'][newPs] += offset
        # Apply labels to boundary pores (trim leading 'pores' if present)
        label = apply_label.split('.')[-1]
        label = 'pore.' + label
        logger.debug('The label \''+label+'\' has been applied')
        self[label] = False
        self[label][newPs] = True

    def add_periodic_connections(self, pores1, pores2, apply_label='periodic'):
        r"""
        Accepts two sets of pores and connects them with new throats.  The
        connections are determined by pairing each pore in ``pores1`` with its
        nearest pore in ``pores2``.  For cubic Networks this will create
        pairings with pores directly across the domain from each other,
        assuming the input pores are 2D co-planar sets of pores.

        Parameters
        ----------
        pores_1 and pores_2 : array_like
            Lists of pores on the opposing faces which are to be linked to
            create periodicity.

        apply_label = string
            The label to apply to the newly created throats.  The default is
            'periodic'.

        Notes
        -----
        This method will raise an exception if the input pores do not create
        fully unique pairs.  Specifically, the length of pore_1 and pores_2
        must be the same AND each pore in pores_1 must pair up with one and
        only one pore in pores_2, and vice versa.  If these conditions are
        not met then periodicity cannot be acheived, and an exception is
        raised.

        """
        logger.debug('Creating periodic pores')
        if sp.shape(pores1)[0] != sp.shape(pores2)[0]:
            raise Exception('Unequal length inputs, periodicity not possible')
        p1 = self['pore.coords'][pores1]
        p2 = self['pore.coords'][pores2]
        dist_mat = sptl.distance_matrix(p1, p2)
        dist_min = sp.amin(dist_mat, axis=1, keepdims=True)
        [a, b] = sp.where(dist_mat == dist_min)
        pairs = sp.vstack([pores1[a], pores2[b]]).T
        # Confirm that each pore in each list is only paired up once
        temp_1 = sp.unique(pairs[:, 0])
        if sp.shape(temp_1) < sp.shape(pores1):
            raise Exception('Non-unique pairs found, periodicity not met')
        temp_2 = sp.unique(pairs[:, 1])
        if sp.shape(temp_2) < sp.shape(pores2):
            raise Exception('Non-unique pairs found, periodicity not met')
        # Add throats to the network for the periodic connections
        self.extend(throat_conns=pairs, labels=apply_label)
        # Create a list which pores are connected which
        self['pore.periodic_neighbor'] = sp.nan
        self['pore.periodic_neighbor'][pairs[:, 0]] = pairs[:, 1]
        self['pore.periodic_neighbor'][pairs[:, 1]] = pairs[:, 0]
        logger.info('Periodic boundary pores added successfully')

    def asarray(self, values):
        r"""
        Retreive values as a rectangular array, rather than the OpenPNM list format

        Parameters
        ----------
        values : array_like
            The values from the network (in a list) to insert into the array

        Notes
        -----
        This method can break on networks that have had boundaries added.  It
        will usually work IF the list of values came only from 'internal' pores.
        """
        if sp.shape(values)[0] > self.num_pores('internal'):
            raise Exception('The received values are bigger than the original'
                            'network')
        Ps = sp.array(self['pore.index'][self.pores('internal')], dtype=int)
        arr = sp.ones(self._shape)*sp.nan
        ind = sp.unravel_index(Ps, self._shape)
        arr[ind[0], ind[1], ind[2]] = values
        return arr

    def fromarray(self, array, propname):
        r"""
        Apply data to the network based on a rectangular array filled with
        values.  Each array location corresponds to a pore in the network.

        Parameters
        ----------
        array : array_like
            The rectangular array containing the values to be added to the
            network. This array must be the same shape as the original network.

        propname : string
            The name of the pore property being added.
        """
        array = sp.atleast_3d(array)
        if sp.shape(array) != self._shape:
            raise Exception('The received array does not match the original network')
        temp = array.flatten()
        Ps = sp.array(self['pore.index'][self.pores('internal')], dtype=int)
        propname = 'pore.' + propname.split('.')[-1]
        self[propname] = sp.nan
        self[propname][self.pores('internal')] = temp[Ps]

    def domain_length(self, face_1, face_2):
        r"""
        Calculate the distance between two faces of the network

        Parameters
        ----------
        face_1 and face_2 : array_like
            Lists of pores belonging to opposite faces of the network

        Returns
        -------
        The length of the domain in the specified direction

        Notes
        -----
        - Does not yet check if input faces are perpendicular to each other
        """
        # Ensure given points are coplanar before proceeding
        if misc.iscoplanar(self['pore.coords'][face_1]) and \
                misc.iscoplanar(self['pore.coords'][face_2]):
            # Find distance between given faces
            x = self['pore.coords'][face_1]
            y = self['pore.coords'][face_2]
            Ds = misc.dist(x, y)
            L = sp.median(sp.amin(Ds, axis=0))
        else:
            logger.warning('The supplied pores are not coplanar. Length will be \
                            approximate.')
            f1 = self['pore.coords'][face_1]
            f2 = self['pore.coords'][face_2]
            distavg = [0, 0, 0]
            distavg[0] = sp.absolute(sp.average(f1[:, 0]) - sp.average(f2[:, 0]))
            distavg[1] = sp.absolute(sp.average(f1[:, 1]) - sp.average(f2[:, 1]))
            distavg[2] = sp.absolute(sp.average(f1[:, 2]) - sp.average(f2[:, 2]))
            L = max(distavg)
        return L

    def domain_area(self, face):
        r"""
        Calculate the area of a given network face

        Parameters
        ----------
        face : array_like
            List of pores of pore defining the face of interest

        Returns
        -------
        The area of the specified face
        """
        coords = self['pore.coords'][face]
        rads = self['pore.diameter'][face]/2.
        # Calculate the area of the 3 principle faces of the bounding cuboid
        dx = max(coords[:, 0]+rads) - min(coords[:, 0]-rads)
        dy = max(coords[:, 1]+rads) - min(coords[:, 1]-rads)
        dz = max(coords[:, 2]+rads) - min(coords[:, 2]-rads)
        yz = dy*dz  # x normal
        xz = dx*dz  # y normal
        xy = dx*dy  # z normal
        # Find the directions parallel to the plane
        directions = sp.where([yz, xz, xy] != max([yz, xz, xy]))[0]
        try:
            # Use the whole network to do the area calculation
            coords = self['pore.coords']
            rads = self['pore.diameter']/2.
            d0 = (max(coords[:, directions[0]]+rads) -
                  min(coords[:, directions[0]]-rads))
            d1 = (max(coords[:, directions[1]]+rads) -
                  min(coords[:, directions[1]]-rads))
            A = d0*d1
        except:
            # If that fails, use the max face area of the bounding cuboid
            A = max([yz, xz, xy])
        if not misc.iscoplanar(self['pore.coords'][face]):
            logger.warning('The supplied pores are not coplanar. Area will be'
                           'approximate')
            pass
        return A

    def subdivide(self, pores=[], shape=[], labels=[]):
        tools.subdivide(network=self, pores=pores, shape=shape, labels=labels)
    subdivide.__doc__ = tools.subdivide.__doc__
