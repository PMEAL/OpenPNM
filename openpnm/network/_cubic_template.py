import logging
import numpy as np
from openpnm.network import Cubic
from openpnm import topotools


logger = logging.getLogger(__name__)
__all__ = ['CubicTemplate']


class CubicTemplate(Cubic):
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

    Examples
    --------
    .. plot::

       import openpnm as op
       import matplotlib.pyplot as plt

       im = op.topotools.template_cylinder_annulus(10, 15, 10)
       pn = op.network.CubicTemplate(template=im)

       # It can be plotted for quick visualization using
       fig, ax = plt.subplots(figsize=(5, 5))
       op.topotools.plot_connections(network=pn, ax=ax, linewidth=0.5)

       plt.show()

    For larger networks and more control over presentation use `Paraview
    <http://www.paraview.org>`_.

    """

    def __init__(self, template, spacing=[1, 1, 1], **kwargs):
        template = np.atleast_3d(template)
        super().__init__(shape=template.shape, spacing=spacing, **kwargs)
        coords = np.unravel_index(range(template.size), template.shape)
        self['pore.template_coords'] = np.vstack(coords).T
        self['pore.template_indices'] = self.Ps
        topotools.trim(network=self, pores=template.flatten() == 0)
        # Add "internal_surface" label to "fake" surface pores!
        ndims = topotools.dimensionality(self).sum()
        max_neighbors = 6 if ndims == 3 else 4
        num_neighbors = np.diff(self.get_adjacency_matrix(fmt="csr").indptr)
        mask_surface = self["pore.surface"]
        mask_internal_surface = (num_neighbors < max_neighbors) & ~mask_surface
        self.set_label("pore.internal_surface", pores=mask_internal_surface)
