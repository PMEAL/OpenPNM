# -*- coding: utf-8 -*-
"""
===============================================================================
CubicDual: Generate a cubic lattice with an interpentrating dual network
===============================================================================

"""
import scipy as sp
from openpnm.network import GenericNetwork, Cubic
from openpnm import topotools
from openpnm.core import logging, Workspace
logger = logging.getLogger(__name__)
ws = Workspace()


class CubicDual(GenericNetwork):
    r"""
    Generates a cubic network of the specified size and shape, then inserts
    a second cubic network into the center of its lattice cells (N-1 pores
    in all directions).  These two networks are further connected by throats
    enabling material to be exchanged between them.

    This network is essentially a *'bcc'* lattice, *except* that the seconary
    network has pores on each face of the domain, which breaks the body-centric
    arranagement.  This allows boundary conditions to be applied to the
    seconary network for transport simuations.

    Parameters
    ----------
    shape : list of ints
        The size and shape of the primary cubic network in terms of the
        number of pores in each direction.  Secondary nodes will be added at
        centers of each unit cell.

    spacing : list of floats
        The distance between pores of the primary network in each of the
        principal directions

    label_1 : string
        The label to apply to the primary cubic lattices, which defaults to
        'primary'

    label_2 : string
        The label to apply to the secondary cubic lattices, which defaults to
        'seconary'

    project : OpenPNM Project object (optional)
        If not provided one will be generated and the network will be assigned
        to it.  It can be retrieved from ``net.project``.

    name : string
        A unique name for the network

    See Also
    --------
    Bravais

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.CubicDual(shape=[3, 3, 3])
    >>> pn.num_pores('pore.primary')  # Normal cubic network is present
    27
    >>> pn.Np  # But more pores are present from seconary network
    59

    And it can be plotted for quick visualization using:

    >>> fig = op.topotools.plot_connections(network=pn,
    ...                                     throats=pn.throats('primary'),
    ...                                     color='b')
    >>> fig = op.topotools.plot_connections(network=pn,
    ...                                     throats=pn.throats('secondary'),
    ...                                     color='r')
    >>> fig = op.topotools.plot_coordinates(network=pn, c='r', s=75, fig=fig)

    .. image:: /../docs/static/images/cubic_dual_network.png
        :align: center

    For larger networks and more control over presentation use `Paraview
    <http://www.paraview.org>`_.

    """
    def __init__(self, shape, spacing=1, label_1='primary',
                 label_2='secondary', **kwargs):
        super().__init__(**kwargs)
        spacing = sp.array(spacing)
        shape = sp.array(shape)
        # Deal with non-3D shape arguments
        shape = sp.pad(shape, [0, 3-shape.size], mode='constant',
                       constant_values=1)
        net = Cubic(shape=shape, spacing=[1, 1, 1])
        net['throat.'+label_1] = True
        net['pore.'+label_1] = True
        single_dim = shape == 1
        shape[single_dim] = 2
        dual = Cubic(shape=shape-1, spacing=[1, 1, 1])
        faces = [['front', 'back'], ['left', 'right'], ['top', 'bottom']]
        faces = [faces[i] for i in sp.where(~single_dim)[0]]
        faces = sp.array(faces).flatten().tolist()
        dual.add_boundary_pores(faces)
        # Add secondary network name as a label
        dual['pore.'+label_2] = True
        dual['throat.'+label_2] = True
        # Shift coordinates prior to stitching
        dual['pore.coords'] += 0.5*(~single_dim)
        topotools.stitch(net, dual, P_network=net.Ps, P_donor=dual.Ps,
                         len_max=1)
        net['throat.interconnect'] = net['throat.stitched']
        del net['throat.stitched']
        net['pore.coords'] *= spacing
        # Clean-up labels
        net['pore.surface'] = False
        net['throat.surface'] = False
        for face in faces:
            # Remove face label from secondary network since it's internal now
            Ps = net.pores(labels=[face, label_2], mode='intersection')
            net['pore.'+face][Ps] = False
            Ps = net.pores(labels=[face+'_boundary'])
            net['pore.'+face][Ps] = True
            Ps = net.pores(face)
            net['pore.surface'][Ps] = True
            Ts = net.find_neighbor_throats(pores=Ps, mode='intersection')
            net['throat.surface'][Ts] = True
            net['throat.'+face] = net.tomask(throats=Ts)
        [net.pop(item) for item in net.labels() if 'boundary' in item]
        # Label non-surface pores and throats as internal
        net['pore.internal'] = ~net['pore.surface']
        Ts = net.find_neighbor_throats(pores=net['pore.internal'])
        net['throat.internal'] = False
        net['throat.internal'][Ts] = True
        # Transfer all dictionary items from 'net' to 'self'
        [self.update({item: net[item]}) for item in net]
        ws.close_project(net.project)
