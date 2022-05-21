import numpy as np
from openpnm.core import LabelMixin, ModelsDict, ModelWrapper, ParserMixin
from openpnm.utils import Workspace, SettingsAttr, PrintableList, PrintableDict
from openpnm.utils import parse_mode
from copy import deepcopy
import inspect
import uuid
import numpy.lib.recfunctions as rf


ws = Workspace()


__all__ = [
    'Base2',
    'Domain',
]


class BaseSettings:
    r"""
    The default settings to use on instance of Base

    Parameters
    ----------
    uuid : str
        A universally unique identifier for the object to keep things straight

    """


class Base2(dict):

    def __init__(self, network=None, settings=None, name='obj'):
        super().__init__()
        # Add settings attribute
        self._settings = SettingsAttr(BaseSettings)
        self.settings._update(settings)
        self.settings['uuid'] = str(uuid.uuid4())
        # Add parameters attr
        self._params = PrintableDict(key="Parameters", value="Value")
        # Associate with project
        if network is None:
            project = ws.new_project()
        else:
            project = network.project
        project.extend(self)
        self.name = name

    def __eq__(self, other):
        return hex(id(self)) == hex(id(other))

    def __repr__(self):
        module = self.__module__
        module = ".".join([x for x in module.split(".") if not x.startswith("_")])
        cname = self.__class__.__name__
        return f'<{module}.{cname} at {hex(id(self))}>'

    def __setitem__(self, key, value):
        if value is None:
            return

        # Intercept parameters
        if key.startswith('param'):
            self._params[key] = value
            return

        # Intercept @ symbol
        if '@' in key:
            element, prop = key.split('@')[0].split('.', 1)
            domain = key.split('@')[1].split('.')[-1]
            locs = super().__getitem__(element+'.'+domain)
            try:
                vals = self[element+'.'+prop]
                vals[locs] = value
                self[element+'.'+prop] = vals
            except KeyError:
                value = np.array(value)
                if value.dtype == bool:
                    temp = np.zeros([self._count(element), *value.shape[1:]],
                                    dtype=bool)
                else:
                    temp = np.zeros([self._count(element), *value.shape[1:]],
                                    dtype=float)*np.nan
                self.__setitem__(element+'.'+prop, temp)
                self[element+'.'+prop][locs] = value
            return

        element, prop = key.split('.', 1)
        # Catch dictionaries and break them up
        if isinstance(value, dict):
            for k, v in value.items():
                self[key+'.'+k] = v
            return
        # Catch dictionaries and convert to struct arrays and back
        # if isinstance(value, dict):
        #     s = self._dict_to_struct(d=value, element=element)
        #     d = self._struct_to_dict(s=s)
        #     for k, v in d.items():
        #         self[element+'.'+prop+'.'+k] = v
        #     return

        # Enfore correct dict naming
        if element not in ['pore', 'throat']:
            raise Exception('All keys must start with either pore, or throat')
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
            try:
                return self._params[key]
            except KeyError:
                return self.network._params[key]

        # If key starts with conduit, then call the get_conduit_data method
        # to build an Nt-by-3 array of pore1-throat-pore2 values
        if key.startswith('conduit'):
            if '@' in key:
                raise Exception('@domain syntax does not work with conduit prefix')
            return self.get_conduit_data(propname=key.split('.', 1)[1])

        # If key contains an @ symbol then return a subset of values at the
        # requested locations, by recursively calling __getitem__
        if '@' in key:
            element, prop = key.split('@')[0].split('.', 1)
            domain = key.split('@')[1].split('.')[-1]
            locs = self[element+'.'+domain]
            vals = self[element+'.'+prop]
            return vals[locs]

        # This allows for lookup of all data that 'ends with' a certain
        # string. Is probably gimmicky and should be deleted.  Was mostly
        # meant for exploring the idea of putting phase values in the
        # network dict, like pn['pore.temperature.air'], but we're not doing
        # that now.
        if '*' in key:
            d = {}
            key = key[1:]
            for k, v in self.items():
                if k.endswith(key):
                    d[k[:-len(key)-1]] = super().__getitem__(k)
            return d

        try:
            return super().__getitem__(key)
        except KeyError:
            # If key is object's name or all, return ones
            if key.split('.', 1)[-1] in [self.name, 'all']:
                element, prop = key.split('.', 1)
                vals = np.ones(self._count(element), dtype=bool)
                self[element+'.all'] = vals
                return vals
            else:
                vals = {}  # Gather any arrays into a dict
                for k in self.keys():
                    if k.startswith(key+'.'):
                        vals.update({k.replace(key+'.', ''): self[k]})
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
                super().__delitem__(key+'.'+item)

    def clear(self, mode=None):
        if mode is None:
            super().clear()
        else:
            mode = parse_mode(obj=self, mode=mode, allowed=['props',
                                                            'labels',
                                                            'models'])
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
            mode = parse_mode(obj=self, mode=mode, allowed=['props',
                                                            'labels',
                                                            'models',
                                                            'constants'])
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
        if self.Np is not None:
            self['pore.'+name] = np.ones([self.Np, ], dtype=bool)
        if self.Nt is not None:
            self['throat.'+name] = np.ones([self.Nt, ], dtype=bool)

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
        # if (self._settings_docs is None) and (settings.__doc__ is not None):
        #     self._settings_docs = settings.__doc__

    def _get_settings(self):
        if self._settings is None:
            self._settings = SettingsAttr()
        # if self._settings_docs is not None:
        #     self._settings.__dict__['__doc__'] = self._settings_docs
        return self._settings

    def _del_settings(self):
        self._settings = None

    settings = property(fget=_get_settings, fset=_set_settings, fdel=_del_settings)

    @property
    def network(self):
        return self.project.network

    def params(self):
        return self._params

    # TODO: Delete this once codes stops asking for it
    @property
    def _domain(self):
        return self

    def _count(self, element):
        for k, v in self.items():
            if k.startswith(element):
                return v.shape[0]

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
            if v.dtype != bool:
                if k.split('.', 1)[0] in element:
                    props.append(k)
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

    def get_conduit_data(self, propname, mode='mean'):
        poreprop = 'pore.' + propname.split('.', 1)[-1]
        throatprop = 'throat.' + propname.split('.', 1)[-1]
        conns = self.network.conns
        try:
            T = self[throatprop]
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

    def _dict_to_struct(self, d, element):
        for k, v in d.items():
            self[element+'._.'+k] = v
        # Do it 2nd loop so that any @ domains are all written before popping
        d2 = {}
        for k, v in d.items():
            # Since key can only be popped once, but may be in items > once
            temp = self.pop(element+'._.'+k.split('@')[0], None)
            if temp is not None:
                d2.update({k.split('@')[0]: temp})
        struct = rf.unstructured_to_structured(np.vstack(list(d2.values())).T,
                                               names=list(d2.keys()))
        return struct

    def _struct_to_dict(self, s):
        d = {}
        for key in s.dtype.names:
            d[key] = s[key]
        return d

    def __str__(self):
        module = self.__module__
        module = ".".join([x for x in module.split(".") if not x.startswith("_")])
        cname = self.__class__.__name__
        horizontal_rule = '―' * 78
        lines = [horizontal_rule]
        lines.append(f"{module}.{cname} : {self.name}")
        lines.append(horizontal_rule)
        lines.append("{0:<5s} {1:<45s} {2:<10s}".format('#',
                                                        'Properties',
                                                        'Valid Values'))
        fmt = "{0:<5d} {1:<45s} {2:>5d} / {3:<5d}"
        lines.append(horizontal_rule)
        props = self.props()
        props.sort()
        for i, item in enumerate(props):
            prop = item
            required = self._count(item.split('.')[0])
            if len(prop) > 35:  # Trim overly long prop names
                prop = prop[0:32] + '...'
            if self[item].dtype == object:  # Print objects differently
                invalid = [i for i in self[item] if i is None]
                defined = np.size(self[item]) - len(invalid)
                lines.append(fmt.format(i + 1, prop, defined, required))
            elif '._' not in prop:
                try:  # normal ndarray
                    a = np.isnan(self[item])
                    defined = np.shape(self[item])[0] \
                        - a.sum(axis=0, keepdims=(a.ndim-1) == 0)[0]
                except TypeError:  # numpy struct array
                    k = self[item].dtype.names
                    required = self[item].shape[0]*len(k)
                    defined = sum([sum(~np.isnan(self[item][i])) for i in k])
                lines.append(fmt.format(i + 1, prop, defined, required))
        lines.append(horizontal_rule)
        # lines = lines.rpartition('\n')[0]
        # lines = lines.append(self._params.__str__())
        return '\n'.join(lines)


class ModelMixin2:

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.models = ModelsDict()

    def add_model(self, propname, model, domain='all', regen_mode='normal',
                  **kwargs):
        if '@' in propname:
            propname, domain = propname.split('@')
        else:
            domain = domain.split('.')[-1]

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

    def add_model_collection(self, models, domain='all'):
        models = deepcopy(models)
        for k, v in models.items():
            _ = v.pop('regen_mode', None)
            model = v.pop('model')
            self.add_model(propname=k, model=model, domain=domain,
                           regen_mode='deferred', **v)

    def regenerate_models(self, propnames=None, exclude=[]):
        if isinstance(propnames, list) and len(propnames) == 0:
            return  # Short-circuit function and return
        elif isinstance(propnames, str):  # Convert string to list if necessary
            propnames = [propnames]
        elif propnames is None:  # If no props given, then regenerate them all
            propnames = self.models.dependency_list()
        # Remove any that are specifically excluded
        propnames = [i for i in propnames if i not in exclude]
        # Re-order given propnames according to dependency tree
        # all_models = self.models.dependency_list()
        all_models = self.models.keys()
        propnames = [i for i in all_models if i in propnames]
        # Now run each on in sequence
        for item in propnames:
            try:
                self.run_model(item)
            except KeyError:
                pass

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
            domain = domain.split('.')[-1]
            element, prop = propname.split('@')[0].split('.', 1)
            propname = element+'.'+prop
            mod_dict = self.models[propname+'@'+domain]
            # Collect kwargs
            kwargs = {'target': self, 'domain': element+'.'+domain}
            for item in mod_dict.keys():
                if item not in ['model', 'regen_mode']:
                    kwargs[item] = mod_dict[item]
            # Deal with models that don't have domain argument yet
            if 'domain' not in inspect.getfullargspec(mod_dict['model']).args:
                _ = kwargs.pop('domain', None)
                vals = mod_dict['model'](**kwargs)
                if isinstance(vals, dict):  # Handle models that return a dict
                    for k, v in vals.items():
                        v = np.atleast_1d(v)
                        if v.shape[0] == 1:  # Returned item was a scalar
                            v = np.tile(v, self._count(element))
                        vals[k] = v[self[element+'.'+domain]]
                elif isinstance(vals, (int, float)):  # Handle models that return a float
                    vals = np.atleast_1d(vals)
                else:  # Index into full domain result for use below
                    vals = vals[self[element+'.'+domain]]
            else:  # Model that accepts domain arg
                vals = mod_dict['model'](**kwargs)
            # Finally add model results to self
            if isinstance(vals, np.ndarray):  # If model returns single array
                if propname not in self.keys():
                    # Create empty array if not found
                    self[propname] = np.nan*np.ones([self._count(element),
                                                     *vals.shape[1:]])
                self[propname][self[element+'.'+domain]] = vals
            elif isinstance(vals, dict):  # If model returns a dict of arrays
                for k, v in vals.items():
                    if propname + '.' + k not in self.keys():
                        # Create empty array if not found
                        self[propname + '.' + k] = \
                            np.nan*np.ones([self._count(element),
                                            *v.shape[1:]])
                    self[propname + '.' + k][self[element+'.'+domain]] = v


class Domain(ParserMixin, LabelMixin, ModelMixin2, Base2):
    ...


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


if __name__ == '__main__':
    import openpnm as op
    import pytest

    # %%
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

    # %% Run some basic tests
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



















