import numpy as np
import pytest
import uuid


class Base:
    r"""
    This class defines all the organizational details like name and domain.
    Everything in OpenPNM will descend from this, and have additional
    functionality added.

    """

    def __init__(self, name=None):
        self.name = name

    def _set_name(self, name):
        if name is None:
            name = str(uuid.uuid4())[:8]
        self._name = name

    def _get_name(self):
        return self._name

    name = property(fget=_get_name, fset=_set_name)

    def _get_domain(self):
        try:
            return self._domain
        except AttributeError:
            return self

    def _set_domain(self, domain):
        self._domain = domain

    domain = property(fset=_set_domain, fget=_get_domain)

    @property
    def network(self): ...

    @property
    def project(self): ...


class NdDataFrame(dict):

    def __init__(self, d={}):
        for k, v in d.items():
            self[k] = v

    def __setitem__(self, key, value):
        N = self.N
        value = np.atleast_1d(value)
        if value.shape[0] == 1:  # Extend scalar to full length
            self._allocate_array(key, shape=(N,), dtype=value.dtype)
            self[key].fill(value[0])
        elif (value.shape[0] == N) or (N == ()):
            super().__setitem__(key, value)
        else:
            raise ValueError(f'Array being assigned to {key} is wrong length')

    def _allocate_array(self, key, dtype=float, shape=[]):
        r"""
        Creates an empty array of the correct size and given type
        """
        N = self.N
        shape = list(shape)
        if N == ():
            raise Exception("NdDataFrame size has not been set yet")

        if shape == []:
            shape = [N, ]

        if shape[0] == N:
            arr = np.zeros(shape, dtype=dtype)
        elif shape[0] < N:
            dtype = bool if dtype == bool else float
            shape[0] = N
            arr = np.zeros(shape, dtype=dtype)
        else:
            raise ValueError(f'Requested size of {key} is too large for NdDataFrame')

        if dtype == float:
            arr.fill(np.nan)

        self[key] = arr

    @property
    def N(self):
        try:
            arr = self[list(self.keys())[0]]
            N = arr.shape[0]
        except:
            N = ()
        return N

    def drop(self, locs):
        mask = np.ones(self.N, dtype=bool)
        mask[locs] = False
        for item in self.keys():
            super().__setitem__(item, self[item][mask])


class DataMixin:
    def __str__(self):
        dom = self.domain
        s = "{"
        for i, key in enumerate(dom.keys()):
            s += " "*(i > 0) + "\'" + key + "\'': " + str(self[key]) + ",\n"
        s += "}"
        return s

    def count(self, key):
        mask = self.locs(key, asmask=True)
        return mask.sum()

    def locs(self, key, asmask=False):
        try:
            element, prop = key.split('.', 1)
        except ValueError:
            element = key
        try:
            locs = self.domain[element + '.' + self.name]
        except KeyError:
            N = self.domain._data[element].N
            locs = np.ones(N, dtype=bool)
        if not asmask:
            locs = np.where(locs)[0]
        return locs

    def set_values(self, key, values, indices=None, relative=True):
        element, prop = key.split('.', 1)
        try:
            arr = self.domain[key]
        except KeyError:
            df = self.domain._data[element]
            temp = np.atleast_1d(values)
            shape = list(temp.shape)
            shape[0] = self.count(key)
            dtype = temp.dtype
            df._allocate_array(key=prop, shape=shape, dtype=dtype)
            arr = self.domain[key]
        locs = self.locs(key)
        if indices is None:
            indices = np.arange(self.count(key))
        if relative:
            arr[locs[indices]] = values
        else:
            arr[indices] = values
        self.domain[key] = arr

    def to_mask(self): ...

    def to_indices(self): ...

    def interpolate_data(self): ...

    def get_conduit_data(self): ...

    def _parse_indices(self): ...

    def _parse_element(self): ...

    def _parse_labels(self): ...

    def _parse_mode(self): ...

    def _parse_prop(self): ...


class LabelMixin:

    def pores(self, label, relative=True):
        return self._get_locs(element='pore', label=label,
                              relative=relative)

    def throats(self, label, relative=True):
        return self._get_locs(element='throat', label=label,
                              relative=relative)

    def _get_locs(self, element, label, relative):
        if label.startswith(element):
            element, label = label.split('.', 1)
        key = element + '.' + label
        if relative:
            locs = np.where(self.domain[key][self.locs(key)])[0]
        else:
            locs = np.where(self.domain[key])[0]
        return locs

    def props(self, element=['pore', 'throat']):
        if type(element) == str:
            element = [element]
        keys = self.domain.keys()
        props = []
        for item in keys:
            if (item.split('.', 1)[0] in element) and (self[item].dtype != bool):
                props.append(item)
        return props

    def labels(self, element=['pore', 'throat']):
        if type(element) == str:
            element = [element]
        keys = self.domain.keys()
        props = []
        for item in keys:
            if (item.split('.', 1)[0] in element) and (self[item].dtype == bool):
                props.append(item)
        return props

    def set_label(self): ...

    def filter_by_label(self): ...

    def num_pores(self): ...

    def num_throats(self): ...


class FullDomain(Base, DataMixin, LabelMixin):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._data = {'pore': NdDataFrame(),
                      'throat': NdDataFrame(),
                      'param': {}}

    def __setitem__(self, key, value):
        element, prop = key.split('.')
        df = self._data[element]  # Get correct data frame
        df[prop] = value

    def __getitem__(self, key):
        element, prop = key.split('.')
        df = self._data[element]  # Get correct data frame
        return df[prop]

    def keys(self):
        keys = []
        for element in self._data.keys():
            for item in self._data[element].keys():
                keys.append(element + '.' + item)
        return keys


class SubDomain(Base, DataMixin, LabelMixin):

    def __init__(self, domain, name, pores=None, throats=None):
        self.domain = domain
        self.name = name
        if pores is not None:
            domain['pore.' + name] = False
            domain['pore.' + name][pores] = True
        if throats is not None:
            domain['throat.' + name] = False
            domain['throat.' + name][throats] = True

    def __setitem__(self, key, value):
        element, prop = key.split('.', 1)
        mask = self.domain[element + '.' + self.name]
        N = self.count(key)
        if np.isscalar(value):
            shape = [self.count(key), ]
        else:
            shape = np.shape(value)
            if shape[0] != N:
                raise ValueError(f'Size of {key} is wrong size for SubDomain')
        try:
            self.domain[key][mask] = value
        except KeyError:
            dtype = bool if np.array(value).dtype == bool else float
            self.domain._data[element]._allocate_array(key=prop,
                                                       shape=shape,
                                                       dtype=dtype)
            self.domain[key][mask] = value

    def __getitem__(self, key):
        element, prop = key.split('.', 1)
        mask = self.domain[element + '.' + self.name]
        val = self.domain[key][mask]
        if len(val) == 0:
            raise KeyError(f'{key}')
        return val

    def set_locations(self, indices, element, mode):
        if mode == 'add':
            self.domain[element + '.' + self.name][indices] = True
        elif mode == 'drop':
            self.domain[element + '.' + self.name][indices] = False
        elif mode == 'clear':
            self.domain[element + '.' + self.name] = False
        elif mode in ['replace', 'overwrite', 'new']:
            self.domain[element + '.' + self.name] = False
            self.domain[element + '.' + self.name][indices] = True
        elif mode == 'purge':
            del self.domain[element + '.' + self.name]


if __name__ == '__main__':
    # Let's test the functionality
    pn = FullDomain()
    pn['pore.blah'] = np.ones(4, dtype=bool)
    geo = SubDomain(domain=pn, name='geom1', pores=[0, 2])
    geo2 = SubDomain(domain=pn, name='geom2', pores=[1, 3])
    # Assign list of values to geo and ensure it is on pn
    geo['pore.test'] = [True, False]
    # Ensure 'test' is on pn, with False in non geo locations
    assert np.all(pn['pore.test'] == [True, False, False, False])
    # Assign scalar and ensure it broadcasts to full length on geo
    geo['pore.bar'] = 2
    assert np.all(geo['pore.bar'] == [2.0, 2.0])
    # And it broadcasts to full length on pn, including 0's
    assert np.isnan(pn['pore.bar']).sum() == 2
    # Assign a boolean to geo and ensure it broadcast's and remains bool
    geo['pore.bool'] = True
    assert geo['pore.bool'].sum() == 2
    assert geo['pore.bool'].dtype == bool
    assert pn['pore.bool'].sum() == 2
    assert pn['pore.bool'].dtype == bool
    # Assign a list of values to pn
    pn['pore.blah'] = [3, 2, 2, 3]
    # Ensure it an be accessed from geo
    assert np.all(geo['pore.blah'] == [3, 2])
    # Assign a scalar to pn and ensure it broadcasts
    pn['pore.test2'] = 1.0
    assert np.all(pn['pore.test2'] == [1.0, 1.0, 1.0, 1.0])
    # Ensure it can be fetched by geo
    assert np.all(geo['pore.test2'] == [1.0, 1.0])
    # Ensure a long list is rejected by geo
    with pytest.raises(ValueError): geo['pore.reject'] = [1, 2, 3]
    with pytest.raises(ValueError): geo['pore.reject'] = [1, 2, 3, 4]
    with pytest.raises(ValueError): geo['pore.reject'] = [1, 2, 3, 4, 5]
    # Ensure key 'reject' wasn't created by accident
    assert 'pore.reject' not in pn.keys()
    # Ensure a short list is rejected by pn
    with pytest.raises(ValueError): pn['pore.reject'] = [1, 2]
    # Ensure key 'reject' wasn't created by accident
    assert 'pore.reject' not in pn.keys()
    # Ensure pn can't accept an array that's too long
    with pytest.raises(ValueError): pn['pore.reject'] = [1, 2, 3, 4, 5]
    # Ensure pn can't accept an array that's too short
    with pytest.raises(ValueError): pn['pore.reject'] = [1, 2, 3]

    print(pn)
    print(geo)
