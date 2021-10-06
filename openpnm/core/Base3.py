import numpy as np
from pprint import pprint
import pytest
import uuid

class Base:
    r"""
    This class defines all the organizational details like name and domain.
    Everything in OpenPNM will descend from this, and have additional
    functionality added.

    """
    def __new__(cls, *args, **kwargs):
        self.settings = {}

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
        self._domian = domain

    domain = property(fset=_set_domain, fget=_get_domain)


class Full(Base, dict, Domain):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        N = self.count(key)
        if np.isscalar(value):
            value = np.ones(N, type(value))*value
        value = np.atleast_1d(value)
        assert value.shape[0] == N, "Received array is wrong length"
        super().__setitem__(key, value)

    def _allocate_array(self, key, value):
        r"""
        Creates an empty array of the correct size and type

        This function is called by the Sub class when attempting to add data
        """
        N = self.count(key)
        value = np.atleast_1d(value)
        dtype = bool if value.dtype == bool else float
        if value.ndim == 1:
            arr = np.zeros(N, dtype=dtype)
        else:
            arr = np.zeros([N, value.shape[1]], dtype=dtype)
        const = {float: np.nan, int: 0, bool: False}
        self[key] = arr*const.get(dtype, None)


class Sub(Base, Data):

    def __init__(self, domain, name=None, **locs):
        super().__init__(name=name)
        self._domain = domain
        for k, v in locs.items():
            mask = np.zeros(domain.count(k), dtype=bool)
            mask[v] = True
            self.domain[k + '.' + name] = mask

    def __setitem__(self, key, value):
        if not np.isscalar(value):
            value = np.atleast_1d(value)
            assert value.shape[0] == self.count(key), "Received array is wrong length"
        try:
            locs = self.locs(key)
            self.domain[key][locs] = value
        except KeyError:
            self.domain._allocate_array(key, value)
            self.__setitem__(key, value)

    def __getitem__(self, key):
        val = self.domain[key][self.locs(key)]
        if len(val) == 0:
            raise KeyError(f'{key}')
        return val


if __name__ == '__main__':
    # Let's test the functionality
    pn = Domain(props={'pore': 4}, name='bob')
    geo = Sub(domain=pn, name='foo', pore=[0, 2])
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
    with pytest.raises(AssertionError): geo['pore.reject'] = [1, 2, 3]
    with pytest.raises(AssertionError): geo['pore.reject'] = [1, 2, 3, 4]
    with pytest.raises(AssertionError): geo['pore.reject'] = [1, 2, 3, 4, 5]
    # Ensure key 'reject' wasn't created by accident
    assert 'pore.reject' not in pn.keys()
    # Ensure a short list is rejected by pn
    with pytest.raises(AssertionError): pn['pore.reject'] = [1, 2]
    # Ensure key 'reject' wasn't created by accident
    assert 'pore.reject' not in pn.keys()
    # Ensure pn can't accept an array that's too long
    with pytest.raises(AssertionError): pn['pore.reject'] = [1, 2, 3, 4, 5]
    # Ensure pn can't accept an array that's too short
    with pytest.raises(AssertionError): pn['pore.5reject'] = [1, 2, 3]


    print(pn)
    print(geo)




















