import numpy as np
from openpnm.network import GenericNetwork
from openpnm import topotools
from openpnm.network.generators import cubic
from openpnm.utils import logging

logger = logging.getLogger(__name__)


class Cubic(GenericNetwork):
    r"""
    Simple cubic lattice with connectivity from 6 to 26

    Parameters
    ----------
    shape : array_like
        The [Nx, Ny, Nz] size of the network in terms of the number of pores in
        each direction
    spacing : array_like, optional
        The spacing between pore centers in each direction. If not given, then
        [1, 1, 1] is assumed.
    connectivity : int, optional
        The number of connections to neighboring pores.  Connections are made
        symmetrically to any combination of face, edge, or corners neighbors.
        The default is 6 to create a simple cubic structure, but options are:

        - 6: Faces only
        - 14: Faces and Corners
        - 18: Faces and Edges
        - 20: Edges and Corners
        - 26: Faces, Edges and Corners

        For a more random distribution of connectivity, use a high
        ``connectivity`` (i.e. 26) and then delete a fraction of the throats
        using ``openpnm.topotools.reduce_coordination``.

    name : string
        An optional name for the object to help identify it.  If not given,
        one will be generated.
    project : OpenPNM Project object, optional
        Each OpenPNM object must be part of a *Project*.  If none is supplied
        then one will be created and this Network will be automatically
        assigned to it.  To create a *Project* use ``openpnm.Project()``.

    """

    def __init__(self, shape, spacing=[1, 1, 1], connectivity=6,
                 name=None, project=None, **kwargs):

        super().__init__(name=name, project=project, **kwargs)

        temp = cubic(shape=shape,
                     spacing=spacing,
                     connectivity=connectivity)
        # fix temp by adding pore and throat prefix
        self.update(temp)
        Np = self['pore.coords'].shape[0]
        Nt = self['throat.conns'].shape[0]
        self['pore.all'] = np.ones((Np, ), dtype=bool)
        self['throat.all'] = np.ones((Nt, ), dtype=bool)
        topotools.label_faces(network=self)
