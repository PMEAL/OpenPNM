import logging
import numpy as np
from openpnm.network import Network
from openpnm._skgraph.generators import cubic_template
from openpnm._skgraph.queries import find_coordination
from openpnm._skgraph.tools import dimensionality, find_surface_nodes


__all__ = ['CubicTemplate']


class CubicTemplate(Network):
    r"""
    Simple cubic lattice with arbitrary domain shape specified by a
    template image

    The class creates a standard Cubic network the same shape as the
    provided image, then trims pores from the network that are not in the
    mask.

    Parameters
    ----------
    template : array_like
        The array (image) describing the desired shape of the domain. All
        locations in the image that are marked as ``True`` are kept while
        the rest of trimmed to yeild the shape.
    spacing : array_like, optional
        The spacing between pore centers in each direction. If not given,
        then [1, 1, 1] is assumed.
    name : str
        An optional name for the object to help identify it. If not given,
        one will be generated.

    Notes
    -----
    The other arguments are the same as ``Cubic`` except that ``shape`` is
    inferred from the ``template`` image.

    """

    def __init__(self, template, spacing=[1, 1, 1], label_surface_pores=False, **kwargs):
        super().__init__(**kwargs)
        template = np.atleast_3d(template)
        net = cubic_template(template=template,
                             spacing=spacing,
                             node_prefix='pore',
                             edge_prefix='throat')
        self.update(net)
        if not label_surface_pores:
            return
        self['pore.surface'] = find_surface_nodes(self)
        ndims = dimensionality(self).sum()
        max_neighbors = 6 if ndims == 3 else 4
        num_neighbors = find_coordination(self, nodes=self.Ps)
        mask_internal_surface = (num_neighbors < max_neighbors) & ~self["pore.surface"]
        self.set_label("pore.internal_surface", pores=mask_internal_surface)
