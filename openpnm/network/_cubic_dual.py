import logging
import numpy as np
from openpnm.network import GenericNetwork, Bravais
from openpnm import topotools
from openpnm.utils import Workspace, Docorator


docstr = Docorator()
logger = logging.getLogger(__name__)
ws = Workspace()
__all__ = ['CubicDual']


@docstr.dedent
class CubicDual(GenericNetwork):
    r"""
    Body centered cubic lattice plus face centered nodes on the surfaces

    This network is essentially a *'bcc'* lattice, *except* that the
    seconary network (body-centered pores) has pores on each face of the
    domain, which breaks the body-centric arranagement. This allows
    boundary conditions to be applied to the seconary network for
    transport simuations.

    Parameters
    ----------
    shape : list[int]
        The size and shape of the primary cubic network in terms of the
        number of pores in each direction. Secondary nodes will be added
        at centers of each unit cell.
    spacing : list[float]
        The distance between pores of the primary network in each of the
        principal directions
    %(GenericNetwork.parameters)s

    See Also
    --------
    Bravais

    """
    def __init__(self, shape, spacing=1, **kwargs):
        super().__init__(**kwargs)
        shape = np.array(shape)
        net = Bravais(shape=shape+1, mode='bcc')
        Ps = net['pore.surface'] * net['pore.corner_sites']

    def add_boundary_pores(self, labels=['top', 'bottom', 'front', 'back',
                                         'left', 'right'], spacing=None):
        r"""
        Add boundary pores to the specified faces of the network

        Pores are offset from the faces by 1/2 of the given ``spacing``,
        such that they lie directly on the boundaries.

        Parameters
        ----------
        labels : str or list[str]
            The labels indicating the pores defining each face where
            boundary pores are to be added (e.g. 'left' or
            ['left', 'right'])
        spacing : scalar or array_like
            The spacing of the network (e.g. [1, 1, 1]).  This must be
            given since it can be quite difficult to infer from the
            network, for instance if boundary pores have already added to
            other faces.

        """
        spacing = np.array(spacing)
        if spacing.size == 1:
            spacing = np.ones(3)*spacing
        for item in labels:
            Ps = self.pores(item)
            coords = np.absolute(self['pore.coords'][Ps])
            axis = np.count_nonzero(np.diff(coords, axis=0), axis=0) == 0
            offset = np.array(axis, dtype=int)*spacing/2
            if np.amin(coords) == np.amin(coords[:, np.where(axis)[0]]):
                offset = -1*offset
            topotools.add_boundary_pores(network=self, pores=Ps, offset=offset,
                                         apply_label=item + '_boundary')
