from openpnm.pnmlib.core import set_data, get_data, get_prop_data, flatten_dict
from openpnm.pnmlib.inspect import tree, get_params
import numpy as np
import pytest


class GenericTest:

    hr = '―' * 78

    def __init__(self):
        print(self.hr)

    def setup_class(self):
        pass

    def run_all(self):
        self.setup_class()
        for item in self.__dir__():
            if item.startswith('test'):
                print(f"Running test: {item}")
                self.__getattribute__(item)()


class GetDataTests(GenericTest):

    def test_get_data(self):
        arr = np.ones(10, dtype=bool)
        d = {'network': {'pore.all': arr}}
        arr2 = get_data(d['network'], 'pore.all')
        assert arr is arr2

    def test_get_data_with_group_in_key(self):
        arr = np.ones(10, dtype=bool)
        d = {'network/pore.all': arr}
        arr2 = get_data(d, 'network/pore.all')
        assert arr is arr2

    def test_get_data_numeric_keys_are_returned(self):
        d = {'network/pore.all': np.ones(10, dtype=bool)}
        arr = get_data(d, 1.1)
        assert arr == 1.1

    def test_get_data_wildcard_for_prop(self):
        d = {'pore.all': np.ones(10, dtype=bool),
             'pore.test': np.ones(10)}
        d2 = get_data(d, key='pore.*')
        assert 'pore.all' in d2.keys()
        assert 'pore.test' in d2.keys()

    def test_get_data_wildcard_for_group(self):
        arr = np.ones(10, dtype=bool)
        d = {'phase1/pore.all': arr,
             'phase2/pore.all': arr,
             'phase3/pore.all': arr,
             }
        d2 = get_data(d, key='*/pore.all')
        assert d2.keys() == d.keys()

    def test_get_data_wildcard_for_group_and_prop(self):
        arr = np.ones(10, dtype=bool)
        d = {'phase/pore.all': arr,
             'phase/phase2/pore.all': arr,
             'phase/phase3/pore.all': arr,
             }
        d2 = get_data(d, key='*/pore.*')
        assert d.keys() == d2.keys()

    def test_get_data_wildcard_for_element(self):
        arr = np.ones(10, dtype=bool)
        d = {'phase.pore.all': arr,
             'phase/phase2/pore.all': arr,
             'phase/phase3/pore.all': arr,
             }
        d2 = get_data(d, key='*.all')
        assert d.keys() == d2.keys()

    def test_get_data_wildcard_for_element_and_subgroup(self):
        arr = np.ones(10, dtype=bool)
        d = {'phase/pore.all': arr,
             'phase/phase2/pore.all': arr,
             'phase/phase3/pore.all': arr,
             }
        d2 = get_data(d, key='*/*/*.all')
        assert len(d2.keys()) < len(d.keys())

    def test_get_data_wildcard_for_element_and_prop(self):
        arr = np.ones(10, dtype=bool)
        d = {'pore.all': arr,
             'pore.test': arr,
             'throat.all': arr,
             'throat.test': arr,
             'param.test': 2.2,
             }
        d2 = get_data(d, key='*.*')
        assert len(d2.keys()) == 5
        assert 'pore.all' in d2.keys()
        assert 'pore.test' in d2.keys()
        assert 'throat.all' in d2.keys()
        assert 'throat.test' in d2.keys()
        assert 'param.test' in d2.keys()

    def test_get_data_using_slash_to_access_nested_groups(self):
        arr = np.ones(10, dtype=bool)
        vals = np.ones(10, dtype=int)
        d = {'phase/pore.all': arr,
             'phase/pore.test1': vals*1,
             'phase/pore.all': arr,
             'phase/pore.test2': vals*2,
             'phase/phase2/pore.all': arr,
             'phase/phase2/pore.test3': vals*3,
             }
        arr3 = get_data(d, key='phase/phase2/pore.test3')
        assert arr3.sum() == 30
        with pytest.raises(KeyError):
            _ = get_data(d, key='phase/phase3/pore.test2')

    def test_get_data_domain_notation(self):
        arr = np.ones(10, dtype=bool)
        vals = np.ones(10, dtype=int)
        left = np.zeros_like(arr)
        left[:4] = True
        right = np.zeros_like(arr)
        right[-3:] = True
        d = {'phase/pore.all': arr,
             'phase/pore.test1': vals*1,
             'phase/pore.left': left,
             'phase/pore.right': right,
             'phase/phase2/pore.all': arr,
             'phase/phase2/pore.test2': vals*2,
             'phase/phase2/pore.left': left,
             'phase/phase2/pore.right': right,
             'phase/phase2/phase3/pore.all': arr,
             'phase/phase2/phase3/pore.test3': vals*3,
             }
        arr2 = get_data(d, key='phase/pore.test1@left')
        assert len(arr2) == sum(left)
        arr3 = get_data(d, key='phase/phase2/pore.test2@right')
        assert len(arr3) == sum(right)
        arr3 = get_data(d, key='phase/phase2/pore.test2@right')
        assert len(arr3) == sum(right)
        with pytest.raises(KeyError):
            _ = get_data(d, key='phase/phase2/pore.test2@top')
        with pytest.raises(KeyError):
            _ = get_data(d, key='phase/phase3/pore.test2@left')

    def test_get_data_params(self):
        arr = np.ones(10, dtype=bool)
        d = {'pore.all': arr,
             'param.MW': 29.1}
        val = get_data(d, 'param.MW')
        assert val == 29.1

    def test_get_data_attr(self):
        arr = np.ones(10, dtype=bool)
        d = {'network/pore.all': arr}
        a = {'attr.ID': 123, 'attr.parent': 321}
        set_data(d, 'network', a)
        d2 = get_data(d, 'network/attr.*')
        assert 'network/attr.ID' in d2.keys()
        assert 'network/attr.parent' in d2.keys()
        d2 = get_data(d, 'network/attr.*')
        assert 'network/attr.ID' in d2.keys()
        assert 'network/attr.parent' in d2.keys()

    def test_get_data_attr_as_prefix(self):
        arr = np.ones(10, dtype=bool)
        d = {'network/pore.all': arr}
        a = {'attr.ID': 123, 'attr.parent': 321}
        set_data(d, 'network', a)
        assert get_data(d, 'network/attr.ID') == 123

    def test_get_data_attr_as_prefix_with_wildcard(self):
        arr = np.ones(10, dtype=bool)
        d = {'pore.all': arr}
        a = {'attr.ID': 123, 'attr.parent': 321}
        set_data(d, '', a)
        attrs = get_data(d, 'attr.*')
        assert attrs == a

    def test_get_data_fetch_attr_from_multiple_groups_with_wildcard(self):
        arr = np.ones(10, dtype=bool)
        d = {}
        set_data(d, 'network/pore.all', arr)
        set_data(d, 'phase/pore.all', arr)
        a = {'attr.ID': 456}
        set_data(d, 'phase', a)
        a = {'attr.ID': 123, 'attr.parent': 321}
        set_data(d, 'network', a)
        attrs = get_data(d, '*/attr.*')
        assert len(attrs) == 3

    def test_get_data_return_target_on_just_wildcard(self):
        arr = np.ones(10, dtype=bool)
        d = {'network/pore.all': arr}
        a = {'attr.ID': 123, 'attr.parent': 321}
        set_data(d, 'network', a)
        d3 = get_data(d, '*')
        assert d3 == d

    # def test_get_prop_data_wildcard_element_and_prop(self):
    #     d = {'pore.all': np.ones(10, dtype=bool),
    #          'pore.test': np.ones(10, dtype=int),
    #          'throat.all': np.ones(20, dtype=bool),
    #          'throat.test': np.ones(20, dtype=float),
    #          'param.test': 2.2,
    #          }
    #     d2 = get_prop_data(d, key='*.*')
    #     assert len(d2.keys()) == 2
    #     assert 'pore.test' in d2.keys()
    #     assert 'throat.test' in d2.keys()
    #     d3 = get_prop_data(d, key='*.*')
    #     assert len(d3.keys()) == 2
    #     assert 'pore.test' in d3.keys()
    #     assert 'throat.test' in d3.keys()

    def test_get_data_components_with_user_defined_delimeter(self):
        phase = {}
        set_data(phase, 'pore.all', np.ones(10, dtype=bool))
        set_data(phase, 'attr.components', ['#o2', '#n2'])
        set_data(phase, 'pore.temperature' + phase['attr.components'][0], 298)
        set_data(phase, 'pore.temperature#n2', 298)
        T = get_data(phase, 'pore.temperature#*')
        assert T['pore.temperature#n2'][0] == 298.0
        assert T['pore.temperature#n2'].size == 10
        assert T['pore.temperature#o2'][0] == 298.0
        assert T['pore.temperature#o2'].size == 10

    def test_get_data_with_pipe(self):
        vals = np.ones(4)
        d = {'phase/pore.all': vals,
             'phase/pore.test1': vals*1,
             'phase/phase2/pore.test2': vals*2,
             'phase/phase2/phase3/pore.test3': vals*3}
        d2 = get_data(d, 'phase|*')
        assert list(d2.keys()) == ['pore.all',
                                   'pore.test1',
                                   'phase2/pore.test2',
                                   'phase2/phase3/pore.test3']
        d2 = get_data(d, 'phase/phase*|*')
        assert list(d2.keys()) == ['pore.test2', 'phase3/pore.test3']


if __name__ == '__main__':

    t = GetDataTests()
    t.run_all()
