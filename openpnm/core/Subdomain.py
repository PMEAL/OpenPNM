from openpnm.core import Base
import numpy as np


class Subdomain(Base):
    r"""
    This subclass of the Base class provides the ability assign the object
    to specific locations (pores and throats) in the domain.  This class
    is subclassed by GenericGeometry and GenericPhysics.

    Notes
    -----
    The following table list the two methods added to Base by this subclass.

    +---------------------+---------------------------------------------------+
    | Methods             | Description                                       |
    +=====================+===================================================+
    | ``add_locations``   | Specified which pores and throats the object      |
    |                     | should be assigned to                             |
    +---------------------+---------------------------------------------------+
    | ``drop_locations``  | Removes the object from the specified pores and   |
    |                     | throats                                           |
    +---------------------+---------------------------------------------------+
    | ``_get_boss``       | Geomtry and Physics objects are subservient to    |
    |                     | Network and Phase objects, respectively, so this  |
    |                     | method returns a handle to the appropriate *boss* |
    +---------------------+---------------------------------------------------+

    Also listed above is a hidden method that might be useful.  The act of
    assign a Subdomain object to a subset of pores or throats basically amounts
    to creating a list in the *boss* object with the Subdomain's name, like
    ``'pore.geo_1'``, with True values where ``geo_1`` applies and False
    elsewhere.  Changing the locations of objects is just a matter of changing
    the locations of the True's and False's.

    The Project object has two methods, ``check_geometry_health`` and
    ``check_physics_health`` that look to make sure all locations are assigned
    to one and only one Geometry and/or Physics.

    """

    def __getitem__(self, key):
        element = key.split('.')[0]
        # Find boss object (either phase or network)
        boss = self.project.find_full_domain(self)
        # Get values if present, or regenerate them
        vals = self.get(key)
        # If still not found, check with boss object
        if vals is None:
            inds = boss._get_indices(element=element, labels=self.name)
            vals = boss[key][inds]
        return vals

    def __setitem__(self, key, value):
        # If value is a dict, skip all this.  The super class will parse
        # the dict individually, at which point the below is called.
        if self.project and not hasattr(value, 'keys'):
            proj = self.project
            boss = proj.find_full_domain(self)
            keys = boss.keys(mode='all', deep=True)
            # Prevent 'pore.foo' on subdomain when already present on boss
            if key in set(boss.keys()).difference(set(self.keys())):
                hit = [i for i in keys if i.startswith(key)][0]
                raise Exception('Cannot create ' + key + ' when ' +
                                hit + ' is already defined')
        super().__setitem__(key, value)

    def add_locations(self, pores=[], throats=[]):
        r"""
        Adds associations between an object and its boss object at the
        given pore and/or throat locations.

        Parameters
        ----------
        pores and throats : array_like
            The pore and/or throat locations for which the association should
            be added.  These indices are for the full domain.

        Notes
        -----
        For *Physics* objects, the boss is the *Phase* with which it was
        assigned, while for *Geometry* objects the boss is the *Network*.

        """
        boss = self.project.find_full_domain(self)
        pores = boss._parse_indices(pores)
        throats = boss._parse_indices(throats)
        if len(pores) > 0:
            self._set_locations(element='pore', indices=pores, mode='add')
        if len(throats) > 0:
            self._set_locations(element='throat', indices=throats, mode='add')

    def drop_locations(self, pores=[], throats=[], complete=False):
        r"""
        Removes association between an objectx and its boss object at the
        given pore and/or throat locations.

        Parameters
        ----------
        pores and throats : array_like
            The pore and/or throat locations from which the association should
            be removed.  These indices refer to the full domain.

        complete : boolean (default is ``False``)
            If ``True`` then *all* pore and throat associations are removed
            along with any trace that the objects were associated.

        Notes
        -----
        For *Physics* objects, the boss is the *Phase* with which it was
        assigned, while for *Geometry* objects the boss is the *Network*.

        """
        boss = self.project.find_full_domain(self)
        pores = boss._parse_indices(pores)
        throats = boss._parse_indices(throats)
        if complete:
            pores = boss.pores(self.name)
            throats = boss.throats(self.name)
        if len(pores) > 0:
            self._set_locations(element='pore', indices=pores, mode='drop',
                                complete=complete)
        if len(throats) > 0:
            self._set_locations(element='throat', indices=throats, mode='drop',
                                complete=complete)

    def _set_locations(self, element, indices, mode, complete=False):
        r"""
        This private method is called by ``set_locations`` and
        ``remove_locations`` as needed.

        """
        boss = self.project.find_full_domain(self)
        element = self._parse_element(element=element, single=True)

        # Make sure label array exists in boss
        if (element+'.'+self.name) not in boss.keys():
            boss[element+'.'+self.name] = False

        # Check to ensure indices aren't already assigned
        if mode == 'add':
            if self._isa('geometry'):
                objs = self.project.geometries().keys()
            else:
                objs = self.project.physics().keys()
            for name in objs:
                if element+'.'+name in boss.keys():
                    if np.any(boss[element+'.'+name][indices]):
                        raise Exception('Given indices are already assigned ' +
                                        'to ' + name)

        # Find mask of existing locations (network indexing)
        mask = boss[element+'.'+self.name]
        # Update mask with new locations (either add or remove)
        if mode == 'add':
            mask = mask + boss._tomask(indices=indices, element=element)
        elif mode == 'drop':
            mask = mask * (~boss._tomask(indices=indices, element=element))
        # Change size of all arrays on self
        for item in self.keys(element=element, mode='all'):
            self.update({item: boss[item][mask]})
        # Update label array in network
        boss[element+'.'+self.name] = mask
        # Remove label from boss if ALL locations are removed
        if mode == 'drop':
            if ~np.any(boss[element+'.'+self.name]):
                if complete:
                    del boss[element+'.'+self.name]
                else:
                    boss[element+'.'+self.name] = False
