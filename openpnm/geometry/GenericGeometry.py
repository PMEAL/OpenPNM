import scipy as sp
from openpnm.core import Base, ModelsMixin
from openpnm.core import logging, Workspace
logger = logging.getLogger(__name__)
ws = Workspace()


class GenericGeometry(Base, ModelsMixin):
    r"""
    GenericGeometry - Base class to construct a Geometry object

    Parameters
    ----------
    network : OpenPNM Network Object

    pores and/or throats : array_like
        The list of pores and throats where this physics applies. If either are
        left blank this will apply the Geometry nowhere.

    name : string
        A unique name to apply to the object.  This name will also be used as a
        label to identify where this Geometry applies.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> Ps = pn.pores()  # Get all pores
    >>> Ts = pn.throats()  # Get all throats
    >>> geom = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)
    """

    def __init__(self, network, pores=[], throats=[], name=None):
        super().__init__(name=name, simulation=network.simulation)
        logger.name = self.name

        # Initialize a label dictionary in the associated network
        network['pore.'+self.name] = False
        network['throat.'+self.name] = False
        self['pore.all'] = sp.ones(shape=sp.size(pores), dtype=bool)
        self['throat.all'] = sp.ones(shape=sp.size(throats), dtype=bool)
        network['pore.'+self.name] = False
        network['pore.'+self.name][pores] = True
        network['throat.'+self.name] = False
        network['throat.'+self.name][throats] = True
        network.simulation.add_geometry(self)
