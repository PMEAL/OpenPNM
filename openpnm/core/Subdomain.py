from openpnm.core import Base


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
        if vals is None:
            inds = boss._get_indices(element=element, labels=self.name)
            vals = boss[key][inds]
        return vals

    def add_locations(self, pores=[], throats=[]):
        r"""
        """
        if len(pores) > 0:
            self._set_locations(element='pore', indices=pores, mode='add')
        if len(throats) > 0:
            self._set_locations(element='throat', indices=throats, mode='add')

    def drop_locations(self, pores=[], throats=[]):
        r"""
        """
        if len(pores) > 0:
            self._set_locations(element='pore', indices=pores, mode='drop')
        if len(throats) > 0:
            self._set_locations(element='throat', indices=throats, mode='drop')

    def _set_locations(self, element, indices=[], mode='add'):
        r"""
        """
        if self._isa('physics'):
            boss = self.project.find_phase(self)
        if self._isa('geometry'):
            boss = self.project.network
        element = self._parse_element(element=element, single=True)
        # Use the network's _parse_indices, since indicies could be 'network'
        # length boolean masks
        indices = boss._parse_indices(indices)

        # Add self's label to boss if not present
        if element+'.'+self.name not in boss.keys():
            boss[element+'.'+self.name] = False

        # Find mask of existing locations (network indexing)
        mask = boss[element+'.'+self.name]
        # Update mask with new locations (either add or remove)
        if mode == 'add':
            mask = mask + boss._tomask(indices=indices, element=element)
        elif mode == 'drop':
            mask = mask ^ (boss._tomask(indices=indices, element=element))
        # Change size of all arrays on self
        for item in self.keys(element=element, mode='all'):
            self.update({item: boss[item][mask]})
        # Update label array in network
        boss[element+'.'+self.name] = mask
