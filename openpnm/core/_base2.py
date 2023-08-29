import numpy as np
import logging
import uuid
from copy import deepcopy
from openpnm.core import (
    LabelMixin,
    ParserMixin,
    ModelsMixin2
)
from openpnm.utils import (
    Workspace,
    SettingsAttr,
    PrintableList,
    PrintableDict,
    Docorator,
    get_printable_props,
    get_printable_labels,
)


docstr = Docorator()
logger = logging.getLogger(__name__)
ws = Workspace()


__all__ = [
    'Base2',
    'Domain',
]


@docstr.get_sections(base='BaseSettings', sections=docstr.all_sections)
@docstr.dedent
class BaseSettings:
    r"""
    The default settings to use on instance of Base

    Parameters
    ----------
    uuid : str
        A universally unique identifier for the object to keep things straight

    """
    default_domain = 'domain_1'


@docstr.get_sections(base='Base', sections=['Parameters'])
@docstr.dedent
class Base2(dict):
    r"""
    A subclassed dictionary used for storing data

    Parameters
    ----------
    network : dict
        An OpenPNM Network object, which is a subclass of a dict
    """

    def __new__(cls, *args, **kwargs):
        instance = super(Base2, cls).__new__(cls, *args, **kwargs)
        # It is necessary to set the SettingsAttr here since some classes
        # use it before calling super.__init__()
        instance.settings = SettingsAttr()
        instance.settings['uuid'] = str(uuid.uuid4())
        return instance

    def __init__(self, network=None, project=None, name='obj_?'):
        super().__init__()
        # Add default settings
        self.settings._update(BaseSettings())
        # Add parameters attr
        self._params = PrintableDict(key="Parameters", value="Value")
        # Associate with project
        if (network is None) and (project is None):
            project = ws.new_project()
        elif project is None:
            project = network.project
        project.append(self)
        self.name = name

    def __eq__(self, other):
        return hex(id(self)) == hex(id(other))

    def __repr__(self):  # pragma: no cover
        module = self.__module__
        module = ".".join([x for x in module.split(".") if not x.startswith("_")])
        cname = self.__class__.__name__
        return f'{self.name} : <{module}.{cname} at {hex(id(self))}>'

    def __setitem__(self, key, value):
        if value is None:
            return

        # Intercept parameters
        if key.startswith('param'):
            _, key = key.split('.', 1)
            self._params[key] = value
            return

        if not (key.startswith('pore.') or key.startswith('throat.')):
            raise Exception("All dict names must start with pore, throat, or param")

        # Intercept @ symbol
        if '@' in key:
            element, prop = key.split('@')[0].split('.', 1)
            domain = key.split('@')[1]
            locs = super().__getitem__(f'{element}.{domain}')
            try:
                vals = self[f'{element}.{prop}']
                vals[locs] = value
                self[f'{element}.{prop}'] = vals
            except KeyError:
                value = np.array(value)
                temp = self._initialize_empty_array_like(value, element)
                self.__setitem__(f'{element}.{prop}', temp)
                self[f'{element}.{prop}'][locs] = value
            return

        element, prop = key.split('.', 1)
        # Catch dictionaries and break them up
        if isinstance(value, dict):
            for k, v in value.items():
                self[f'{key}.{k}'] = v
            return

        # Enforce correct dict naming
        if element not in ['pore', 'throat']:
            raise Exception('All keys must start with either pore or throat')
        # Convert value to ndarray
        if not isinstance(value, np.ndarray):
            value = np.array(value, ndmin=1)
        # Skip checks for coords and conns
        if key in ['pore.coords', 'throat.conns']:
            self.update({key: value})
            return
        # Finally write data
        if self._count(element) is None:
            self.update({key: value})  # If length not defined, do it
        elif value.shape[0] == 1:  # If value is scalar
            value = np.ones((self._count(element), ), dtype=value.dtype)*value
            self.update({key: value})
        elif np.shape(value)[0] == self._count(element):
            self.update({key: value})
        else:
            raise Exception('Provided array is wrong length for ' + key)

    def __getitem__(self, key):
        # If key is a just a numerical value, then kick it directly back.
        # This allows one to do either value='pore.blah' or value=1.0 in
        # pore-scale models
        if not isinstance(key, str):
            return key

        if key.startswith('param'):
            _, key = key.split('.', 1)
            try:
                return self._params[key]
            except KeyError:
                return self.network._params[key]

        # If key contains an @ symbol then return a subset of values at the
        # requested locations, by recursively calling __getitem__
        if '@' in key:
            element, prop = key.split('@')[0].split('.', 1)
            domain = key.split('@')[1]
            if f'{element}.{domain}' not in self.keys():
                raise KeyError(key)
            locs = self[f'{element}.{domain}']
            vals = self[f'{element}.{prop}']
            return vals[locs]

        try:
            return super().__getitem__(key)
        except KeyError:
            # If key is object's name or all, return ones
            if key.split('.', 1)[-1] in [self.name, 'all']:
                element, prop = key.split('.', 1)
                vals = np.ones(self._count(element), dtype=bool)
                return vals
            else:
                vals = {}  # Gather any arrays into a dict
                for k in self.keys():
                    if k.startswith(f'{key}.'):
                        vals.update({k.replace(f'{key}.', ''): self[k]})
                if len(vals) > 0:
                    return vals
                else:
                    raise KeyError(key)

    def __delitem__(self, key):
        try:
            super().__delitem__(key)
        except KeyError:
            d = self[key]  # If key is a nested dict, get all values
            for item in d.keys():
                super().__delitem__(f'{key}.{item}')

    def pop(self, *args):
        r"""
        """
        v = super().pop(*args)
        if v is None:
            try:
                d = self[args[0]]
                v = {}
                for item in d.keys():
                    key = f'{args[0]}.{item}'
                    v[key] = super().pop(key)
            except KeyError:
                pass
        return v

    def clear(self, mode=None):
        r"""
        Clears or deletes certain things from object. If no arguments are provided
        it defaults to the normal `dict` behavior.

        Parameters
        ----------
        mode : str
            Controls which things are to be deleted. Options are:

            =========== ============================================================
            `mode`        Description
            =========== ============================================================
            'props'     Deletes all pore and throat properties (i.e numerical data)
                        in the object's dictionary (except 'pore.coords' and
                        'throat.conns' if it is a network object).
            'labels'    Deletes all labels (i.e. boolean data) in the object's
                        dictionary.
            'models'    Delete are pore and throat properties that were produced
                        by a pore-scale model.
            =========== ============================================================

        """
        if mode is None:
            super().clear()
        else:
            if isinstance(mode, str):
                mode = [mode]
            if 'props' in mode:
                for item in self.props():
                    if item not in ['pore.coords', 'throat.conns']:
                        del self[item]
            if 'labels' in mode:
                for item in self.labels():
                    if item not in ['pore.'+self.name, 'throat.'+self.name]:
                        del self[item]
            if 'models' in mode:
                for item in self.models.keys():
                    _ = self.pop(item.split('@')[0], None)

    def keys(self, mode=None):
        r"""
        An overloaded version of ``keys`` that optionally accepts a ``mode``

        Parameters
        ----------
        mode : str
            If given, optionally, it controls which type of keys are returned.
            Options are:

            ========== =======================================================
            mode       description
            ========== =======================================================
            props      Returns only keys that contain numerical arrays
            labels     Returns only keys that contain boolean arrays
            models     Returns only keys that were generated by a pore-scale
                       model
            constants  Returns only keys are were *not* generated by a pore-
                       scale model
            ========== =======================================================

        """
        if mode is None:
            return super().keys()
        else:
            if isinstance(mode, str):
                mode = [mode]
            vals = set()
            if 'props' in mode:
                for item in self.props():
                    vals.add(item)
            if 'labels' in mode:
                for item in self.labels():
                    vals.add(item)
            if 'models' in mode:
                for item in self.models.keys():
                    propname = item.split('@')[0]
                    if propname in self.keys():
                        vals.add(propname)
            if 'constants' in mode:
                vals = vals.union(set(self.props()))
                for item in self.models.keys():
                    propname = item.split('@')[0]
                    if propname in vals:
                        vals.remove(propname)
            return PrintableList(vals)

    def _set_name(self, name, validate=True):
        if not hasattr(self, '_name'):
            self._name = ''
        old_name = self._name
        if name == old_name:
            return
        name = self.project._generate_name(name)
        self._name = name

    def _get_name(self):
        try:
            return self._name
        except AttributeError:
            return None

    name = property(_get_name, _set_name)

    def _get_project(self):
        for proj in list(ws.values()):
            if self in proj:
                return proj

    project = property(fget=_get_project)

    def _set_settings(self, settings):
        self._settings = deepcopy(settings)

    def _get_settings(self):
        if self._settings is None:
            self._settings = SettingsAttr()
        return self._settings

    def _del_settings(self):
        self._settings = None

    settings = property(fget=_get_settings, fset=_set_settings, fdel=_del_settings)

    @property
    def network(self):
        r"""
        Shortcut to retrieve a handle to the network object associated with the
        calling object
        """
        return self.project.network

    @property
    def params(self):
        r"""
        This attribute stores 'scalar' data that can be used by pore-scale models.
        For instance, if a model calls for `temperature` you can specify
        `pore.temperature` if every pore might have a different value, or
        `param.temperature` if a single value prevails everywhere.
        """
        return self._params

    def _count(self, element):
        if element == 'pore':
            try:
                return self['pore.coords'].shape[0]
            except KeyError:
                for k, v in self.items():
                    if k.startswith('pore.'):
                        return v.shape[0]
        elif element == 'throat':
            try:
                return self['throat.conns'].shape[0]
            except KeyError:
                for k, v in self.items():
                    if k.startswith('throat.'):
                        return v.shape[0]

    @property
    def Nt(self):
        r"""
        Shortcut to retrieve the number of throats in the domain
        """
        return self._count('throat')

    @property
    def Np(self):
        r"""
        Shortcut to retrieve the number of pores in the domain
        """
        return self._count('pore')

    @property
    def Ts(self):
        r"""
        Shortcut to retrieve the indices of *all* throats
        """
        return np.arange(self._count('throat'))

    @property
    def Ps(self):
        r"""
        Shortcut to retrieve the indices of *all* pores
        """
        return np.arange(self._count('pore'))

    def _tomask(self, element, indices):
        return self.to_mask(**{element+'s': indices})

    def to_mask(self, pores=None, throats=None):
        r"""
        Generates a boolean mask with `True` values in the given locations

        Parameters
        ----------
        pores : array_like
            The pore indices where `True` values will be placed. If `pores` is
            given the `throats` is ignored.
        throats : array_like
            The throat indices where `True` values will be placed. If `pores` is
            given the `throats` is ignored.

        Returns
        -------
        mask : ndarray, boolean
            A boolean array of length Np is `pores` was given or Nt if
            `throats` was given.

        """
        if pores is not None:
            indices = np.array(pores, ndmin=1)
            N = self.Np
        elif throats is not None:
            indices = np.array(throats, ndmin=1)
            N = self.Nt
        else:
            raise Exception('Must specify either pores or throats')
        mask = np.zeros((N, ), dtype=bool)
        mask[indices] = True
        return mask

    def to_indices(self, mask):
        r"""
        Converts a boolean mask to pore or throat indices

        Parameters
        ----------
        mask : ndarray
            A boolean mask with `True` values indicating either pore or
            throat indices. This array must either be Nt or Np long, otherwise
            an Exception is raised.

        Returns
        -------
        indices : ndarray
            An array containing numerical indices of where `mask` was `True`.

        Notes
        -----
        This function is equivalent to just calling `np.where(mask)[0]` but
        does check to ensure that `mask` is a valid length.
        """
        mask = np.array(mask, dtype=bool)
        if mask.shape[0] not in [self.Np, self.Nt]:
            raise Exception('Mask must be either Nt or Np long')
        return np.where(mask)[0]

    def props(self, element=['pore', 'throat']):
        r"""
        Retrieves a list of keys that contain numerical data (i.e. "properties")

        Parameters
        ----------
        element : str, list of strings
            Indicates whether `'pore'` or `'throat'` properties should be returned.
            The default is `['pore', 'throat']`, so both are returned.

        Returns
        -------
        props : list of strings
            The names of all dictionary keys on the object that contain
            numerical data.
        """
        if element is None:
            element = ['pore', 'throat']
        if isinstance(element, str):
            element = [element]
        props = []
        for k, v in self.items():
            el, prop = k.split('.', 1)
            if (el in element) and (v.dtype != bool) and not prop.startswith('_'):
                props.append(k)
        props = sorted(props)
        props = PrintableList(props)
        return props

    def interpolate_data(self, propname, mode='mean'):
        r"""
        Generates an array of the requested pore/throat data by interpolating
        the neighboring throat/pore data.

        Parameters
        ----------
        propname : str
            The data to be generated.
        mode : str
            Dictate how the interpolation is done. Options are 'mean', 'min',
            and 'max'.

        Returns
        -------
        data : ndarray
            An ndarray containing the interpolated data.  E.g. Requesting
            'throat.temperature' will read the values of 'pore.temperature'
            in each of the neighboring pores and compute the average
            (if `mode='mean'`).
        """
        from openpnm.models.misc import from_neighbor_throats, from_neighbor_pores
        element, prop = propname.split('.', 1)
        if element == 'throat':
            if self['pore.'+prop].dtype == bool:
                raise Exception('The requested datatype is boolean, cannot interpolate')
            values = from_neighbor_pores(self, prop='pore.'+prop, mode=mode)
        elif element == 'pore':
            if self['throat.'+prop].dtype == bool:
                raise Exception('The requested datatype is boolean, cannot interpolate')
            values = from_neighbor_throats(self, prop='throat.'+prop, mode=mode)
        return values

    def get_conduit_data(self, propname):
        r"""
        Fetches an Nt-by-3 array of the requested property

        Parameters
        ----------
        propname : str
            The dictionary key of the property to fetch.

        Returns
        -------
        data : ndarray
            An Nt-by-3 array with each column containing the requrested data
            for pore1, throat, and pore2 respectively.

        """
        poreprop = 'pore.' + propname.split('.', 1)[-1]
        throatprop = 'throat.' + propname.split('.', 1)[-1]
        conns = self.network.conns
        try:
            T = self[throatprop]
            if T.ndim > 1:
                raise Exception(f'{throatprop} must be a single column wide')
        except KeyError:
            T = np.ones([self.Nt, ], dtype=float)*np.nan
        try:
            P1, P2 = self[poreprop][conns.T]
        except KeyError:
            P1 = np.ones([self.Nt, ], dtype=float)*np.nan
            P2 = np.ones([self.Nt, ], dtype=float)*np.nan
        vals = np.vstack((P1, T, P2)).T
        if np.isnan(vals).sum() == vals.size:
            raise KeyError(f'{propname} not found')
        return vals

    def __str__(self):  # pragma: no cover
        hr = '―' * 78
        lines = ''
        lines += '\n' + "═"*78 + '\n' + self.__repr__() + '\n' + hr + '\n'
        lines += get_printable_props(self)
        lines += '\n' + hr + '\n'
        lines += get_printable_labels(self)
        lines += '\n' + hr
        return lines

    def _initialize_empty_array_like(self, value, element):
        element = element.split('.', 1)[0]
        value = np.array(value)
        if value.dtype == bool:
            temp = np.zeros([self._count(element), *value.shape[1:]],
                            dtype=bool)
        else:
            temp = np.zeros([self._count(element), *value.shape[1:]],
                            dtype=float)*np.nan
        return temp


class Domain(ParserMixin, LabelMixin, ModelsMixin2, Base2):
    r"""
    The main class used for Network, Phase and Algorithm objects.

    This class inherits from ``Base2``, but also has several mixins for
    added functionality.

    """
    ...
