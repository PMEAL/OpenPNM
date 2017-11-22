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

    def set_locations(self, pores=[], throats=[], mode='add'):
        r"""
        Assign or unassign a Geometry object to specified locations

        Parameters
        ----------
        pores : array_like
            The pore locations in the Network where this Geometry is to apply

        throats : array_like
            The throat locations in the Network where this Geometry is to apply

        mode : string
            Either 'add' (default) or 'remove' the object from the specified
            locations

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.Np
        125
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,
        ...                                         pores=sp.arange(5, 125),
        ...                                         throats=pn.Ts)
        >>> [geom.Np, geom.Nt]
        [120, 300]
        >>> geom['pore.dummy'] = True
        >>> health = pn.check_geometry_health()
        >>> pores = health['undefined_pores']
        >>> geom.set_locations(pores=pores)
        >>> [geom.Np, geom.Nt]
        [125, 300]
        The label 'pore.dummy' was assigned 'before' these pores were added
        >>> geom.pores(labels='dummy', mode='not')
        array([0, 1, 2, 3, 4])
        >>> geom.set_locations(pores=pores, mode='remove')
        >>> [geom.Np, geom.Nt]
        [120, 300]
        # All pores without 'pore.dummy' label are gone
        >>> geom.num_pores(labels='dummy', mode='not')
        0
        """
        if mode == 'add':
            # Check if any constant values exist on the object
            for item in self.props():
                if (item not in self.models.keys()) or \
                   (self.models[item]['regen_mode'] == 'constant'):
                    raise Exception('Constant properties found on object, ' +
                                    'cannot increase size')
            if sp.size(pores) > 0:
                Tools.SetLocations.add(obj=self, element='pore',
                                       locations=pores)
            if sp.size(throats) > 0:
                Tools.SetLocations.add(obj=self, element='throat',
                                       locations=throats)
        if mode == 'remove':
            if sp.size(pores) > 0:
                Tools.SetLocations.drop(obj=self, element='pore',
                                        locations=pores)
            if sp.size(throats) > 0:
                Tools.SetLocations.drop(obj=self, element='throat',
                                        locations=throats)
        # Finally, regenerate models to correct the length of all arrays
        self.models.regenerate()
