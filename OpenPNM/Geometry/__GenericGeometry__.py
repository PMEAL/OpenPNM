# -*- coding: utf-8 -*-
"""
===============================================================================
GenericGeometry -- Base class to manage pore scale geometry
===============================================================================

"""
import scipy as sp
from OpenPNM.Base import Core
from OpenPNM.Postprocessing import Plots
from OpenPNM.Base import logging, Tools
from OpenPNM.Network import GenericNetwork
logger = logging.getLogger(__name__)


class GenericGeometry(Core):
    r"""
    GenericGeometry - Base class to construct a Geometry object

    Parameters
    ----------
    network : OpenPNM Network Object

    pores and/or throats : array_like
        The list of pores and throats where this physics applies. If either are
        left blank this will apply the Geometry nowhere.  The locations can be
        changed after instantiation using ``set_locations()``.

    name : string
        A unique name to apply to the object.  This name will also be used as a
        label to identify where this geometry applies.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> Ps = pn.pores()  # Get all pores
    >>> Ts = pn.throats()  # Get all throats
    >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,
    ...                                         pores=Ps,
    ...                                         throats=Ts)
    """

    def __init__(self, network=None, pores=[], throats=[], **kwargs):
        super().__init__(**kwargs)
        logger.name = self.name

        if network is None:
            network = GenericNetwork()
        self.network.update({network.name: network})  # Attach network to self
        # Register self with network.geometries
        self._net.geometries.update({self.name: self})

        # Initialize a label dictionary in the associated network
        self._net['pore.'+self.name] = False
        self._net['throat.'+self.name] = False
        try:
            self.set_locations(pores=pores, throats=throats)
        except:
            self.workspace.purge_object(self)
            raise Exception('Provided locations are in use, instantiation cancelled')

    def __getitem__(self, key):
        element = key.split('.')[0]
        # Convert self.name into 'all'
        if key.split('.')[-1] == self.name:
            key = element + '.all'

        if key in list(self.keys()):  # Look for data on self...
            return super(GenericGeometry, self).__getitem__(key)
        if key == 'throat.conns':  # Handle specifically
            [P1, P2] = \
                self._net['throat.conns'][self._net[element+'.'+self.name]].T
            Pmap = sp.zeros((self._net.Np,), dtype=int) - 1
            Pmap[self._net.pores(self.name)] = self.Ps
            conns = sp.array([Pmap[P1], Pmap[P2]]).T
            # Replace -1's with nans
            if sp.any(conns == -1):
                conns = sp.array(conns, dtype=object)
                conns[sp.where(conns == -1)] = sp.nan
            return conns
        else:  # ...Then check Network
            return self._net[key][self._net[element+'.'+self.name]]

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

    def plot_histograms(self,
                        throat_diameter='throat.diameter',
                        pore_diameter='pore.diameter',
                        throat_length='throat.length'):

        Plots.distributions(obj=self,
                            throat_diameter=throat_diameter,
                            pore_diameter=pore_diameter,
                            throat_length=throat_length)

    plot_histograms.__doc__ = Plots.distributions.__doc__
