from openpnm.core import Base
import numpy as np


class Subdomain(Base):
    r"""
    This subclass of the Base class provides the ability assign the object
    to specific locations (pores and throats) in the domain.  This class
    is subclassed by GenericGeometry and GenericPhysics.

    Notes
    -----
    The Project object has two methods, ``check_geometry_health`` and
    ``check_physics_health`` that look to make sure all locations are
    assigned to one and only one Geometry and/or Physics.

    """

    def __getitem__(self, key):
        element = key.split('.')[0]
        # Try to get vals directly first
        vals = self.get(key)
        if vals is None:  # Otherwise invoke search
            # Find boss object (either phase or network)
            boss = self.project.find_full_domain(self)
            inds = boss._get_indices(element=element, labels=self.name)
            try:  # Will invoke interleave data if necessary
                vals = boss[key]  # Will return nested dict if present
                if isinstance(vals, dict):  # Index into each array in nested dict
                    for item in vals:
                        vals[item] = vals[item][inds]
                else:  # Otherwise index into single array
                    vals = vals[inds]
            except KeyError:
                vals = super().__getitem__(key)
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
                raise Exception('Cannot create ' + key + ' when '
                                + hit + ' is already defined')
        super().__setitem__(key, value)

    def _set_locations(self, element, indices, mode):
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
                        raise Exception('Given indices already assigned to ' + name)

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
                boss[element+'.'+self.name] = False
