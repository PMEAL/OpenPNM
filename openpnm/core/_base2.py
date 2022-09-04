import numpy as np
import inspect
import logging
import uuid
from copy import deepcopy
from openpnm.core import (
    LabelMixin,
    ModelsDict,
    ModelWrapper,
    ParserMixin,
)
from openpnm.utils import (
    Workspace,
    SettingsAttr,
    PrintableList,
    PrintableDict,
    Docorator,
    get_printable_props,
    get_printable_labels,
    prettify_logger_message
)


docstr = Docorator()
logger = logging.getLogger(__name__)
ws = Workspace()


__all__ = [
    'Base2',
    'ModelMixin2',
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
    Base class

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

    def __repr__(self):
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
                if value.dtype == bool:
                    temp = np.zeros([self._count(element), *value.shape[1:]],
                                    dtype=bool)
                else:
                    temp = np.zeros([self._count(element), *value.shape[1:]],
                                    dtype=float)*np.nan
                self.__setitem__(f'{element}.{prop}', temp)
                self[f'{element}.{prop}'][locs] = value
            return

        element, prop = key.split('.', 1)
        # Catch dictionaries and break them up
        if isinstance(value, dict):
            for k, v in value.items():
                self[f'{key}.{k}'] = v
            return

        if prop == 'all':
            logger.warning('Cannot use all as a label')

        # Enfore correct dict naming
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

        # If key starts with 'conduit.', then call the get_conduit_data method
        # to build an Nt-by-3 array of pore1-throat-pore2 values
        if key.startswith('conduit'):
            domain = None
            if '@' in key:
                key, domain = key.split('@')
            vals = self.get_conduit_data(propname=key.split('.', 1)[1])
            if domain is not None:
                locs = self['throat.'+domain]
                vals = vals[locs]
            return vals

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

        # This allows for lookup of all data that 'ends with' a certain
        # string, like pn[*diameter'] will return a dict with pore and throat
        # as the dict keys
        # if key.startswith('*'):
        #     d = {}
        #     key = key[1:]  # Remove astrisk
        #     if key.startswith('.'):  # Remove leading dot if *. was given
        #         key = key[1:]
        #     for k, v in self.items():
        #         if k.endswith(key):
        #             d[k[:-len(key)-1]] = super().__getitem__(k)
        #     return d

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
        for proj in ws.values():
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
        return self.project.network

    @property
    def params(self):
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

    # TODO: Delete this once codes stops asking for it
    @property
    def _domain(self):
        return self

    @property
    def Nt(self):
        return self._count('throat')

    @property
    def Np(self):
        return self._count('pore')

    @property
    def Ts(self):
        return np.arange(self._count('throat'))

    @property
    def Ps(self):
        return np.arange(self._count('pore'))

    def _tomask(self, element, indices):
        return self.to_mask(**{element+'s': indices})

    def to_mask(self, pores=None, throats=None):
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
        mask = np.array(mask, dtype=bool)
        return np.where(mask)[0]

    def props(self, element=['pore', 'throat']):
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
        from openpnm.models.misc import from_neighbor_throats, from_neighbor_pores
        element, prop = propname.split('.', 1)
        if element == 'throat':
            if self['pore.'+prop].dtype == bool:
                raise Exception('The requested datatype is boolean, cannot interpolate')
            values = from_neighbor_pores(target=self, prop='pore.'+prop, mode=mode)
        elif element == 'pore':
            if self['throat.'+prop].dtype == bool:
                raise Exception('The requested datatype is boolean, cannot interpolate')
            values = from_neighbor_throats(target=self, prop='throat.'+prop, mode=mode)
        return values

    def get_conduit_data(self, propname):
        r"""
        Fetches an Nt-by-3 array of the request property in the form of
        [pore1, throat, pore2].

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
        # If requested data is already in the form of conduit data, return it
        if propname in self.keys():
            vals = self[propname]
            if (vals.ndim == 2) and (vals.shape[1] == self.Nt):
                return vals
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

    def __str__(self):
        hr = '―' * 78
        lines = ''
        lines += '\n' + "═"*78 + '\n' + self.__repr__() + '\n' + hr + '\n'
        lines += get_printable_props(self)
        lines += '\n' + hr + '\n'
        lines += get_printable_labels(self)
        lines += '\n' + hr
        return lines


class ModelMixin2:

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.models = ModelsDict()

    def add_model(self, propname, model, domain=None, regen_mode='normal',
                  **kwargs):
        if '@' in propname:
            propname, domain = propname.split('@')
        elif domain is None:
            domain = self.settings['default_domain']
        element, prop = propname.split('.', 1)
        if (element + '.' + domain not in self.keys()) and (domain != 'all'):
            try:
                N = self._count(element)
                self[element + '.' + domain] = True
            except:
                logger.warning(f'Could not define {element}.{domain} since '
                               f'number of {element}s is not defined yet')
        domain = domain.split('.', 1)[-1]

        # Add model and regen_mode to kwargs dictionary
        kwargs.update({'model': model, 'regen_mode': regen_mode})
        # Insepct model to extract arguments and default values
        if model.__defaults__:
            vals = list(inspect.getfullargspec(model).defaults)
            keys = inspect.getfullargspec(model).args[-len(vals):]
            for k, v in zip(keys, vals):  # Put defaults into kwargs
                if k not in kwargs:  # Skip if argument was given in kwargs
                    kwargs.update({k: v})
        self.models[propname+'@'+domain] = ModelWrapper(**kwargs)
        if regen_mode != 'deferred':
            self.run_model(propname+'@'+domain)

    def add_model_collection(self, models, regen_mode='deferred', domain=None):
        models = deepcopy(models)
        for k, v in models.items():
            if 'domain' not in v.keys():
                v['domain'] = domain
            if 'regen_mode' not in v.keys():
                v['regen_mode'] = regen_mode
            self.add_model(propname=k, **v)

    def regenerate_models(self, propnames=None, exclude=[]):
        all_models = self.models.dependency_list()
        # Regenerate all properties by default
        if propnames is None:
            propnames = all_models
        else:
            propnames = np.atleast_1d(propnames).tolist()
        # Remove any that are specifically excluded
        propnames = np.setdiff1d(propnames, exclude).tolist()
        # Reorder given propnames according to dependency tree
        tmp = [e.split("@")[0] for e in propnames]
        idx_sorted = [all_models.index(e) for e in tmp]
        propnames = [elem for i, elem in sorted(zip(idx_sorted, propnames))]
        # Now run each on in sequence
        for item in propnames:
            try:
                self.run_model(item)
            except KeyError as e:
                msg = (f"{item} was not run since the following property"
                       f" is missing: {e}")
                logger.warning(prettify_logger_message(msg))
                self.models[item]['regen_mode'] = 'deferred'

    def run_model(self, propname, domain=None):
        if domain is None:
            if '@' in propname:  # Get domain from propname if present
                propname, domain = propname.split('@')
                self.run_model(propname=propname, domain=domain)
            else:  # No domain means run model for ALL domains
                for item in self.models.keys():
                    if item.startswith(propname+'@'):
                        domain = item.split('@')[-1]
                        self.run_model(propname=propname, domain=domain)
        else:  # domain was given explicitly
            domain = domain.split('.', 1)[-1]
            element, prop = propname.split('@')[0].split('.', 1)
            propname = f'{element}.{prop}'
            mod_dict = self.models[propname+'@'+domain]
            # Collect kwargs
            kwargs = {'domain': f'{element}.{domain}'}
            for item in mod_dict.keys():
                if item not in ['model', 'regen_mode']:
                    kwargs[item] = mod_dict[item]
            # Deal with models that don't have domain argument yet
            if 'domain' not in inspect.getfullargspec(mod_dict['model']).args:
                _ = kwargs.pop('domain', None)
                vals = mod_dict['model'](self, **kwargs)
                if isinstance(vals, dict):  # Handle models that return a dict
                    for k, v in vals.items():
                        v = np.atleast_1d(v)
                        if v.shape[0] == 1:  # Returned item was a scalar
                            v = np.tile(v, self._count(element))
                        vals[k] = v[self[f'{element}.{domain}']]
                elif isinstance(vals, (int, float)):  # Handle models that return a float
                    vals = np.atleast_1d(vals)
                else:  # Index into full domain result for use below
                    vals = vals[self[f'{element}.{domain}']]
            else:  # Model that accepts domain arg
                vals = mod_dict['model'](self, **kwargs)
            # Finally add model results to self
            if isinstance(vals, np.ndarray):  # If model returns single array
                if propname not in self.keys():
                    # Create empty array if not found
                    self[propname] = np.nan*np.ones([self._count(element),
                                                     *vals.shape[1:]])
                self[propname][self[f'{element}.{domain}']] = vals
            elif isinstance(vals, dict):  # If model returns a dict of arrays
                for k, v in vals.items():
                    if f'{propname}.{k}' not in self.keys():
                        # Create empty array if not found
                        self[f'{propname}.{k}'] = \
                            np.nan*np.ones([self._count(element), *v.shape[1:]])
                    self[f'{propname}.{k}'][self[f'{element}.{domain}']] = v


class Domain(ParserMixin, LabelMixin, ModelMixin2, Base2):
    ...


if __name__ == '__main__':
    import pytest
    import openpnm as op

    def random_seed(target, domain, seed=None, lim=[0, 1]):
        inds = target[domain]
        np.random.seed(seed)
        seeds = np.random.rand(inds.sum())*(lim[1]-lim[0]) + lim[0]
        return seeds

    def factor(target, prop, f=1):
        vals = target[prop]*f
        return vals

    def dolittle(target, domain):
        N = target[domain].sum()
        d = {}
        d['item1'] = np.ones([N, ])
        d['item2'] = np.ones([N, ])*2
        d['item3'] = np.ones([N, ])*3
        return d

    g = op.network.Cubic(shape=[3, 3, 1])

    g.add_model(propname='pore.seed@left',
                model=random_seed,
                lim=[0.2, 0.4],
                regen_mode='deferred')
    g.add_model(propname='pore.seed@right',
                model=random_seed,
                lim=[0.7, 0.99],
                regen_mode='deferred')
    g.add_model(propname='pore.seedx',
                model=factor,
                prop='pore.seed',
                f=10,
                regen_mode='deferred')
    g.add_model(propname='pore.test',
                model=factor,
                domain='front',
                prop='pore.seed',
                f=100,
                regen_mode='deferred')
    g.add_model(propname='pore.dict2',
                model=dolittle,
                regen_mode='deferred',)
    g.add_model(propname='pore.dict3',
                model=dolittle,
                domain='left',
                regen_mode='deferred',)
    g.add_model(propname='pore.dict3',
                model=dolittle,
                domain='right',
                regen_mode='deferred',)

    ## Run some basic tests

    # Use official args
    g.run_model('pore.seed', domain='pore.left')
    assert np.sum(~np.isnan(g['pore.seed'])) == g['pore.left'].sum()
    # Use partial syntax
    g.run_model('pore.seed', domain='left')
    assert np.sum(~np.isnan(g['pore.seed'])) == g['pore.left'].sum()
    # Use lazy syntax
    g.run_model('pore.seed@right')
    assert np.sum(np.isnan(g['pore.seed'])) == 3
    # Run the pore seed model for all domains at once:
    del g['pore.seed']
    g.run_model('pore.seed')  # This does not work yet
    assert 'pore.seed' in g.keys()
    assert g['pore.seed@left'].shape[0] == 3
    assert g['pore.seed@right'].shape[0] == 3
    assert 'pore.seedx' not in g.keys()
    # Full domain model
    g.run_model('pore.seedx')
    assert 'pore.seedx' in g.keys()
    x = g['pore.seedx']
    assert x[~np.isnan(x)].min() > 2
    # Run the non-domain enabled model on a subdomain
    assert 'pore.test' not in g
    g.run_model('pore.test')
    assert 'pore.test' in g
    assert np.isnan(g['pore.test']).sum() == 7
    del g['pore.test']
    g.run_model('pore.test@front')
    assert np.isnan(g['pore.test']).sum() == 7
    # Fetch data with lazy syntax
    assert g['pore.seed@left'].shape[0] == 3
    # Write data with lazy syntax, ensuring scalar to array conversion
    g['pore.seed@right'] = np.nan
    assert np.sum(~np.isnan(g['pore.seed'])) == g['pore.left'].sum()
    # Write array directly
    g['pore.seed@right'] = np.ones(3)*3
    assert np.sum(np.isnan(g['pore.seed'])) == 3
    # Use labels that were not used by models
    assert g['pore.seed@front'].shape[0] == 3
    # Write a dict
    g['pore.dict'] = {'item1': 1, 'item2': 2.0}
    assert g['pore.dict.item1'].sum() == 9
    assert g['pore.dict.item2'].sum() == 18
    # A dict with domains
    g['pore.dict'] = {'item1@left': 2, 'item2@right': 3.0}
    assert g['pore.dict.item1'].sum() == 12
    assert g['pore.dict.item2'].sum() == 21
    g['pore.dict'] = {'item3@left': 1, 'item3@right': 2.0}
    # Double dots...not sure how these should work
    g['pore.nested.name1'] = 10
    g['pore.nested.name2'] = 20
    assert isinstance(g['pore.nested'], dict)
    assert len(g['pore.nested']) == 2
    with pytest.raises(KeyError):
        g['pore.nested.fail']
    del g['pore.nested.name1']
    assert 'pore.nested.name1' not in g.keys()
    del g['pore.nested']
    assert 'pore.nested.name2' not in g.keys()
    # More fun
    c = g['conduit.seed']
    assert c.shape == (12, 3)
    assert 'throat.seed' not in g
    assert np.isnan(c[:, 1]).sum() == g.Nt
    # Run model that returns a dict to all pores
    g.run_model('pore.dict2')
    assert len(g['pore.dict2']) == 3
    assert g['pore.dict2.item1'].sum() == 9
    assert g['pore.dict2.item2'].sum() == 18
    # Run model that returns a dict to only subdomain pores
    g.run_model('pore.dict3@left')
    assert len(g['pore.dict3']) == 3
    assert np.isnan(g['pore.dict3.item1']).sum() == 6
    assert g['pore.dict3.item1@left'].sum() == 3
    assert np.isnan(g['pore.dict3.item2']).sum() == 6
    assert np.isnan(g['pore.dict3.item3']).sum() == 6
    # Run model on both domains that its assigned
    g.run_model('pore.dict3')
    assert g['pore.dict3.item1@left'].sum() == 3
    assert g['pore.dict3.item1@right'].sum() == 3
