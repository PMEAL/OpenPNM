import openpnm as op
import scipy as sp
import pytest


class BaseTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])

    def teardown_class(self):
        mgr = op.Base.Workspace()
        mgr.clear()

    def test_pores(self):
        Ps = self.net.pores('top')
        assert sp.all(Ps == [2, 5, 8, 11, 14, 17, 20, 23, 26])

    def test_pores_multiple_labels_union(self):
        Ps = self.net.pores(['top', 'front'], mode='union')
        assert sp.all(Ps == [0, 1, 2, 3, 4, 5, 6, 7, 8,
                             11, 14, 17, 20, 23, 26])

    def test_pores_multiple_labels_intersection(self):
        Ps = self.net.pores(['top', 'front'], mode='intersection')
        assert sp.all(Ps == [2, 5, 8])

    def test_parse_indices_boolean(self):
        b = sp.array([True, True, True])
        with pytest.raises(Exception):
            self.net._parse_indices(indices=b)
        b = sp.zeros((self.net.Np,), dtype=bool)
        assert len(self.net._parse_indices(indices=b)) == 0
        b = sp.zeros((self.net.Nt,), dtype=bool)
        b[[0, 1, 2]] = True
        assert sp.shape(self.net._parse_indices(indices=b)) == (3,)

    def test_parse_indices_None(self):
        assert len(self.net._parse_indices(indices=None)) == 0

    def test_parse_indices_string(self):
        with pytest.raises(Exception):
            self.net._parse_indices(indices='abc')

    def test_parse_indices_int(self):
        a = self.net._parse_indices(indices=0)
        assert type(a) == sp.ndarray
        assert sp.all(a == 0)

    def test_parse_indices_list(self):
        a = self.net._parse_indices(indices=[0, 1])
        assert type(a) == sp.ndarray
        assert sp.all(a == [0, 1])

    def test_parse_element_None(self):
        a = self.net._parse_element(element=None)
        assert sorted(a) == ['pore', 'throat']

    def test_parse_element_various_strings(self):
        a = self.net._parse_element(element='pore')
        assert a == ['pore']
        a = self.net._parse_element(element='Pore')
        assert a == ['pore']
        a = self.net._parse_element(element='pores')
        assert a == ['pore']
        a = self.net._parse_element(element='Pores')
        assert a == ['pore']
        a = self.net._parse_element(element='throat')
        assert a == ['throat']
        a = self.net._parse_element(element='Throat')
        assert a == ['throat']
        a = self.net._parse_element(element='throats')
        assert a == ['throat']
        a = self.net._parse_element(element='Throats')
        assert a == ['throat']

    def test_parse_element_bad_string(self):
        with pytest.raises(Exception):
            self.net._parse_element(element='pore2')

    def test_parse_element_duplicate(self):
        a = self.net._parse_element(element=['pore', 'pore'])
        assert a == ['pore']
        a = self.net._parse_element(element=['pore', 'pore'], single=True)
        assert a == 'pore'

    def test_parse_element_single_true(self):
        with pytest.raises(Exception):
            self.net._parse_element(element=['pore', 'throat'], single=True)
        a = self.net._parse_element(element=['pore'], single=True)
        assert a == 'pore'

    def test_parse_element_props(self):
        a = self.net._parse_element(element=['pore.diameter'], single=True)
        assert a == 'pore'

    def test_parse_labels_none(self):
        with pytest.raises(Exception):
            self.net._parse_labels(labels=None, element='pore')

    def test_parse_labels_string(self):
        a = self.net._parse_labels(labels='top', element='pore')
        assert a == ['pore.top']
        a = self.net._parse_labels(labels='internal', element='throat')
        assert a == ['throat.internal']
        a = self.net._parse_labels(labels='pore.top', element='pore')
        assert a == ['pore.top']
        a = self.net._parse_labels(labels='throat.internal', element='throat')
        assert a == ['throat.internal']

    def test_parse_labels_wildcards(self):
        a = self.net._parse_labels(labels='pore.b*', element='pore')
        assert sorted(a) == ['pore.back', 'pore.bottom']
        a = self.net._parse_labels(labels='pore.*ight', element='pore')
        assert sorted(a) == ['pore.right']

    def test_parse_labels_duplicates(self):
        a = self.net._parse_labels(['pore.r*', 'pore.right'], element='pore')
        assert a == ['pore.right']

    def test_parse_mode_string(self):
        a = self.net._parse_mode(mode='union')
        assert a == ['union']

    def test_parse_mode_single(self):
        a = self.net._parse_mode(mode=['union', 'intersection'])
        assert sorted(a) == ['intersection', 'union']
        with pytest.raises(Exception):
            a = self.net._parse_mode(mode=['union1', 'union2'], single=True)
        a = self.net._parse_mode(mode=['union1'], single=True)
        assert a == 'union1'

    def test_parse_mode_allowed(self):
        allowed = ['a', 'b', 'c']
        with pytest.raises(Exception):
            self.net._parse_mode(mode=['a', 'd'], allowed=allowed)

    def test_parse_mode_duplicate(self):
        a = self.net._parse_mode(mode=['union', 'union'])
        assert a == ['union']
        a = self.net._parse_mode(mode=['union', 'union'], single=True)
        assert a == 'union'


if __name__ == '__main__':

    t = BaseTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
