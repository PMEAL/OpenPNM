from openpnm.core import Base
import numpy as np


class Subdomain(Base):

    def __getitem__(self, key):
        element = key.split('.')[0]
        # Find boss object (either phase or network)
        if self._isa('physics'):
            boss = self.project.find_phase(self)
        else:
            boss = self.project.network
        # Deal with a few special key items
        if key.split('.')[-1] == '_id':
            inds = boss._get_indices(element=element, labels=self.name)
            return boss[element+'._id'][inds]
        # Convert self.name into 'all'
        elif key.split('.')[-1] in [self.name]:
            return self[element+'.all']
        # Now get values if present, or regenerate them
        vals = self.get(key)
        # If still not found, check with boss object (interleave data)
        if vals is None:
            inds = boss._get_indices(element=element, labels=self.name)
            vals = boss[key][inds]
        return vals

    def __setitem__(self, key, value):
        if self.project:
            # Find boss object (either phase or network)
            if self._isa('phase'):
                boss = self.project.find_phase(self)
            else:
                boss = self.project.network
            if key in set(boss.keys()).difference({'pore.all', 'throat.all'}):
                raise Exception(key + ' already exists on ' + boss.name)
        super().__setitem__(key, value)

    def add_locations(self, pores=[], throats=[]):
        r"""
        Adds associations between an objectx and its boss object at the
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
        boss = self._get_boss()
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
        boss = self._get_boss()
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
        """
        boss = self._get_boss()
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

    def _get_boss(self):
        if self._isa('physics'):
            boss = self.project.find_phase(self)
        if self._isa('geometry'):
            boss = self.project.network
        return boss
