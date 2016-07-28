# -*- coding: utf-8 -*-
"""
===============================================================================
CubicDUal: Generate lattice-like networks with a
===============================================================================

"""
import scipy as sp
import OpenPNM.Utilities.misc as misc
from OpenPNM.Network import GenericNetwork
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class CubicDual(GenericNetwork):
    r"""
    Generates a cubic network of the specified size and shape.

    Parameters
    ----------
    name : string
        A unique name for the network

    shape : list of ints
        The size and shape of the principal cubic network in terms of the
        number of pores in each direction.  Secondary nodes will be added at
        the center of each unit cell defined by a cube 8 pores in the primary
        network.

    spacing : list of floats
        The distance between pores in each of the principal directions

    label_1 and label_2 : strings
        The labels to apply to the main cubic lattice and the interpenetrating
        cubic lattice (i.e. the dual network).  The defaults are 'corners' and
        'centers', which refers to the position on the unit cell.

    Examples
    --------
    >>>
    """
    def __init__(self, shape=None, spacing=[1, 1, 1], label_1='corners',
                 label_2='centers', **kwargs):
        super().__init__(**kwargs)
        import OpenPNM as op
        spacing = sp.array(spacing)
        shape = sp.array(shape)
        net = op.Network.Cubic(shape=shape, spacing=[1, 1, 1])
        net['throat.'+label_1] = True
        net['pore.'+label_1] = True
        dual = op.Network.Cubic(shape=shape-1)
        dual['pore.'+label_2] = True
        dual['throat.'+label_2] = True
        dual['pore.coords'] += 0.5
        op.Network.tools.stitch(net, dual, P_network=net.Ps,
                                P_donor=dual.Ps, len_max=1)
        net['throat.interconnect'] = net['throat.stitched']
        del net['throat.stitched']
        net['pore.surface'] = False
        net['throat.surface'] = False
        # Create center pores on each surface
        offset = {'back': [0.5, 0, 0], 'front': [-0.5, 0, 0],
                  'right': [0, 0.5, 0], 'left': [0, -0.5, 0],
                  'top': [0, 0, 0.5], 'bottom': [0, 0, -0.5]}
        surface_labels = ['top', 'bottom', 'left', 'right', 'front', 'back']
        for item in surface_labels:
            # Clone 'center' pores and shift to surface
            Ps = net.pores(labels=[item, label_2], mode='intersection')
            net.clone_pores(pores=Ps, apply_label=[label_2, 'surface', item])
            Ps = net.pores(labels=['surface', item], mode='intersection')
            net['pore.coords'][Ps] += offset[item]
            Ts = net.find_neighbor_throats(pores=Ps)
            # Label pores and throats
            net['pore.surface'][Ps] = True
            net['throat.'+label_2][Ts] = True
        Ps = net.pores(labels=surface_labels)
        net['pore.surface'][Ps] = True

        net['pore.coords'] *= spacing
        [self.update({item: net[item]}) for item in net]
        del self.workspace[net.name]

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
