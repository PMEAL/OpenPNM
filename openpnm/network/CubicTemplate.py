import scipy as sp
from openpnm.network import Cubic
from openpnm import topotools
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class CubicTemplate(Cubic):
    r"""
    Simple cubic lattice with arbitrary domain shape specified by a template
    image

    The class creates a standard Cubic network the same shape as the provided
    image, then trims pores from the network that are not in the mask.

    Parameters
    ----------
    template : array_like
        The array (image) describing the desired shape of the domain.  All
        locations in the image that are marked as ``True`` are kept while the
        rest of trimmed to yeild the shape.

    spacing : array_like, optional
        The spacing between pore centers in each direction. If not given, then
        [1, 1, 1] is assumed.

    name : string
        An optional name for the object to help identify it.  If not given,
        one will be generated.

    project : OpenPNM Project object, optional
        Each OpenPNM object must be part of a *Project*.  If none is supplied
        then one will be created and this Network will be automatically
        assigned to it.  To create a *Project* use ``openpnm.Project()``.

    Notes
    -----
    The other arguments are the same as ``Cubic`` except that ``shape`` is
    inferred from the ``template`` image.

    See Also
    --------
    The following methods in ``topotools`` can help generate template images:

    template_cylinder_annulus
    template_sphere_shell

    Examples
    --------
    >>> import openpnm as op
    >>> im = op.topotools.template_cylinder_annulus(15, 10, 5)
    >>> pn = op.network.CubicTemplate(template=im)

    And it can be plotted for quick visualization using:

    >>> fig = op.topotools.plot_connections(network=pn)

    .. image:: /../docs/static/images/cubic_template_network.png
        :align: center

    For larger networks and more control over presentation use `Paraview
    <http://www.paraview.org>`_.

    """
    def __init__(self, template, spacing=[1, 1, 1], **kwargs):

        template = sp.atleast_3d(template)
        if 'shape' in kwargs:
            del kwargs['shape']
            logger.warning('shape argument ignored, inferred from template')
        super().__init__(shape=template.shape, spacing=spacing, **kwargs)

        coords = sp.unravel_index(range(template.size), template.shape)
        self['pore.template_coords'] = sp.vstack(coords).T
        self['pore.template_indices'] = self.Ps
        self['pore.drop'] = template.flatten() == 0
        topotools.trim(network=self, pores=self.pores('drop'))
        del self['pore.drop']
        # remove labels pertaining to surface pores, then redo post-trim
        self.clear(mode='labels')
        self['pore.internal'] = True
        self['throat.internal'] = True
        topotools.find_surface_pores(self)
