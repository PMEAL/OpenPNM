from openpnm.core import Base
import numpy as np


class Subdomain(Base):

    def __getitem__(self, key):
        # Find boss object (either phase or network)
        element = key.split('.')[0]
        if self._isa('phase'):
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

    def add_locations(self, pores=[], throats=[]):
        r"""
        """
        boss = self._get_boss()
        pores = boss._parse_indices(pores)
        throats = boss._parse_indices(throats)
        if len(pores) > 0:
            self._set_locations(element='pore', indices=pores, mode='add')
        if len(throats) > 0:
            self._set_locations(element='throat', indices=throats, mode='add')

    def drop_locations(self, pores=[], throats=[]):
        r"""
        """
        boss = self._get_boss()
        pores = boss._parse_indices(pores)
        throats = boss._parse_indices(throats)
        if len(pores) > 0:
            self._set_locations(element='pore', indices=pores, mode='drop')
        if len(throats) > 0:
            self._set_locations(element='throat', indices=throats, mode='drop')

    def _set_locations(self, element, indices=[], mode='add'):
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
                del boss[element+'.'+self.name]

    def _get_boss(self):
        if self._isa('physics'):
            boss = self.project.find_phase(self)
        if self._isa('geometry'):
            boss = self.project.network
        return boss
