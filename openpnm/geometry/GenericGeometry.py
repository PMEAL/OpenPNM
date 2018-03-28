import scipy as sp
import matplotlib.pyplot as plt
from openpnm.core import Subdomain, ModelsMixin, Workspace, logging
logger = logging.getLogger(__name__)
ws = Workspace()


class GenericGeometry(Subdomain, ModelsMixin):
    r"""
    GenericGeometry - Base class to construct a Geometry object

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network object to which this Geometry applies.

    pores and/or throats : array_like
        The list of pores and throats where this Geometry applies.

    name : string
        A unique name to apply to the object.  This name will also be used as a
        label to identify where this Geometry applies.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> Ps = pn.pores('all')  # Get all pores
    >>> Ts = pn.throats('all')  # Get all throats
    >>> geom = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)
    """

    def __init__(self, network=None, project=None, pores=[], throats=[],
                 settings={}, **kwargs):
        # Define some default settings
        self.settings.update({'prefix': 'geo'})
        # Overwrite with user supplied settings, if any
        self.settings.update(settings)

        # Deal with network or project arguments
        if network is not None:
            project = network.project

        super().__init__(project=project, **kwargs)

        if network is not None:
            network['pore.'+self.name] = False
            network['throat.'+self.name] = False

        self.add_locations(pores=pores, throats=throats)

    def show_hist(self, props=['pore.diameter'], bins=20, fig=None, **kwargs):
        if fig is None:
            fig = plt.figure()
        if type(props) is str:
            props = [props]
        N = len(props)
        if N == 1:
            r = 1
            c = 1
        elif N < 4:
            r = 1
            c = N
        else:
            r = int(sp.ceil(N**0.5))
            c = int(sp.floor(N**0.5))

        for i in range(len(props)):
            plt.subplot(r, c, i+1)
            plt.hist(self[props[i]], bins=bins, **kwargs)

    @property
    def network(self):
        r"""
        A shortcut to get a handle to the associated network
        There can only be one so this works
        """
        return self.project.network
