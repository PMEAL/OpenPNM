from openpnm.core import Base, LabelMixin, ParamMixin
import numpy as np


class Subdomain(Base, LabelMixin):
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
        Element agnostic version of ``set_locations``

        Parameters
        ----------
        element : str
            Whether 'pore' or 'throat' values are being set
        indices : array_like
            The global indices of the pores or throats being set
         mode : str
            Controls how the locations are set.  Options are:

            * 'switch'
                Assigned the object to the specified locations and unassigns
                any other objects from those locations if present
            * 'add'
                Assigns the object to the specified locations. An error is
                riased if the locations are already assigned to another object
            * 'drop'
                Removes the object from the specified locations

        Notes
        -----
        This was originally kept for backward compatibility since
        ``_set_locations`` is called thoughout the code; however, it turned
        out to be necessary to have a generic version that accepted
        ``element``.

        """
        if element.startswith('pore'):
            self.set_locations(pores=indices, mode=mode)
        if element.startswith('throat'):
            self.set_locations(throats=indices, mode=mode)

    def set_locations(self, pores=None, throats=None, mode='switch'):
        r"""
        Assign a Subdomain object to specific pores and/or throats

        Parameters
        ----------
        pores : array_like
            The list of pore indices to which this Subdomain should be
            assigned
        throats : array_like
            The list of throat indices to which this Subdomain should be
            assigned
        mode : str
            Controls how the locations are set.  Options are:

            * 'switch' (default)
                Assigns the object to the specified locations and unassigns
                any other objects from those locations if present
            * 'add'
                Assigns the object to the specified locations. An error is
                raised if the locations are already assigned to another object
            * 'drop'
                Removes the object from the specified locations

        """
        self._parse_mode(mode=mode, allowed=['add', 'drop', 'switch'])
        if (pores is not None) and (throats is not None):
            # If both are sent call method for each then return
            self.set_locations(pores=pores, mode=mode)
            self.set_locations(throats=throats, mode=mode)
            return
        elif pores is not None:
            indices = pores
            element = 'pore'
        elif throats is not None:
            indices = throats
            element = 'throat'

        boss = self._domain

        # Make sure label array exists in boss
        if (element + '.' + self.name) not in boss.keys():
            boss[element + '.' + self.name] = False

        # Find mask of existing locations (global indexing)
        mask = boss[element + '.' + self.name]

        # Update mask with new locations (either add or remove)
        if mode == 'add':
            # Check to ensure indices aren't already assigned
            for i in boss._subdomains:
                if element + '.' + i.name in boss.keys():
                    if np.any(boss[element + '.' + i.name][indices]):
                        raise Exception(f'Indices already assigned to {i.name}',
                                        ', use mode = switch instead')
            # If no exception was raised, generate a mask
            mask = mask + boss._tomask(indices=indices, element=element)
        elif mode == 'drop':
            mask = mask * (~boss._tomask(indices=indices, element=element))
        elif mode == 'switch':
            for i in boss._subdomains:
                i._set_locations(element=element, indices=indices, mode='drop')
            mask = mask + boss._tomask(indices=indices, element=element)

        # Change size of all arrays on self
        for item in self.keys(element=element, mode='all'):
            self.update({item: boss[item][mask]})

        # Finally assign mask to boss
        boss[element + '.' + self.name] = np.copy(mask)

    def to_global(self, pores=None, throats=None):
        r"""
        Convert local indices from a subdomain object to global values

        Parameters
        ----------
        pores, throats : array_like
            List of pore or throat indices to be converted

        Returns
        -------
        indices : ndarray
            An array of location indices
        """
        if pores is not None:
            element = 'pore'
            locs = pores
        elif throats is not None:
            element = 'throat'
            locs = throats
        mask = self._domain[element + '.' + self.name]
        inds = np.where(mask)[0]
        return inds[locs]

    def to_local(self, pores=None, throats=None, missing_vals=-1):
        r"""
        Convert global indices to local values relative to a subdomain object

        Parameters
        ----------
        pores, throats : array_like
            List of pore or throat indices to be converted
        missing_values : scalar
            The value to put into missing locations if global indices are not
            found.  If ``missing_vals`` is ``None``, then any missing values
            are removed from the returned list.

        Returns
        -------
        indices : ndarray
            An array of location indices
        """
        if pores is not None:
            element = 'pore'
            locs = pores
        if throats is not None:
            element = 'throat'
            locs = throats
        mask = np.ones_like(self._domain[element + '.all'],
                            dtype=int)*-1
        inds = np.where(self._domain[element + '.' + self.name])[0]
        mask[inds] = self._get_indices(element)
        vals = mask[locs]
        if missing_vals is None:
            vals = vals[vals >= 0]
        else:
            vals[vals >= 0] = missing_vals
        return vals
