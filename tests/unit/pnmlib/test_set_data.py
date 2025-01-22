from openpnm.pnmlib.core import set_data, get_data
from openpnm.pnmlib.inspect import tree
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


class SetDataTests(GenericTest):

    def test_set_data_with_various_prop_name_and_dict_key_combos(self):
        pn = {}
        # Specify number of pores by assigning a 'pore' array to top level
        set_data(pn, 'pore.all', np.ones(10, dtype=bool))
        # Now scalars can be broadcast to the correct length
        set_data(pn, 'pore.test.one', 1.0)
        assert pn['pore.test.one'].sum() == 10
        # Keys can have extra dots in them
        set_data(pn, 'pore.test.two', 2.0)
        assert pn['pore.test.two'].sum() == 20
        # Dictionaries can be passed as values to write multiple things at once
        set_data(pn, 'network', {'pore.test.zero': 0.0})
        assert pn['network/pore.test.zero'].sum() == 0
        set_data(pn, '', {'pore.test.three': 3.0, 'pore.test.four': 4.0})
        assert pn['pore.test.three'].sum() == 30
        assert pn['pore.test.four'].sum() == 40

    def test_set_data_prop_name_dict_key_edge_cases(self):
        arr = np.ones(10, dtype=int)
        pn = {}
        set_data(pn, 'pore.all', np.ones(10, dtype=bool))

        # Passing an empty dict does nothing
        k = list(pn.keys())
        set_data(pn, 'pore.test', {})
        assert k == list(pn.keys())

        # A new group can be created on the fly, using dicts
        set_data(pn, 'phase1', {'pore.nonsense': arr})
        assert pn['phase1/pore.nonsense'].sum() == 10
        set_data(pn, '', {'phase3/pore.nonsense': arr})
        assert pn['phase3/pore.nonsense'].sum() == 10
        # Scalars are expanded to length of the top level array, even for subgroups
        set_data(pn, 'phase2', {'pore.nonsense': 1})  # No length yet
        assert pn['phase2/pore.nonsense'].sum() == 10

        # Writing must include a valid dict name somewhere
        with pytest.raises(Exception):
            set_data(pn, 'phase2', {'phase3': 1})
        # Replace 1 with a valid dict and voila
        set_data(pn, 'phase2', {'phase3': {'pore.all': arr}})
        assert pn['phase2/phase3/pore.all'].sum() == 10
        set_data(pn, 'phase2', {'phase3/phase4': {'pore.all': arr}})
        assert pn['phase2/phase3/phase4/pore.all'].sum() == 10
        set_data(pn, 'phase2', {'phase3': {'phase4/pore.test': 2*arr}})
        assert pn['phase2/phase3/phase4/pore.test'].sum() == 20
        # Duplicate group names just create an other layer
        set_data(pn, '', {'phase2': {'phase3/pore.test3': 3*arr}})
        assert pn['phase2/phase3/pore.test3'].sum() == 30

    def test_set_data_check_protected_prefixes(self):
        arr = np.ones(10, dtype=bool)
        pn = {}
        set_data(pn, 'pore.all', arr)
        set_data(pn, '', {'pore.test': arr})
        assert pn['pore.test'].sum() == 10
        with pytest.raises(Exception):
            set_data(pn, 'pore/pore.test', arr)
        with pytest.raises(Exception):
            set_data(pn, 'attr', {'attr.test': 'ID'})

    def test_set_data_with_various_group_names_and_dicts(self):
        pn = {}
        set_data(pn, 'phase2', {'pore.one': np.ones(10),
                                'pore.two': 2*np.ones(10)})
        assert pn['phase2/pore.one'].sum() == 10
        assert pn['phase2/pore.two'].sum() == 20
        set_data(pn, '', {'phase2/pore.three': 3*np.ones(10),
                          'phase2/pore.four': 4*np.ones(10)})
        assert pn['phase2/pore.three'].sum() == 30
        assert pn['phase2/pore.four'].sum() == 40
        set_data(pn, '', {'phase2/phase3/pore.five': 5*np.ones(10),
                          'phase2/pore.six': 6*np.ones(10)})
        assert pn['phase2/phase3/pore.five'].sum() == 50
        assert pn['phase2/pore.six'].sum() == 60

    def test_set_data_with_group_name_and_group_in_prefix(self):
        pn = {}
        set_data(pn, 'phase2', {'phase3/pore.five': 5*np.ones(10),
                                'pore.six': 6*np.ones(10)})
        assert pn['phase2/phase3/pore.five'].sum() == 50
        assert pn['phase2/pore.six'].sum() == 60

    def test_set_data_using_slash_in_key(self):
        d = {'pore.all': np.ones(10, dtype=bool)}
        arr = np.ones(10)
        d.update({'pore.test': arr})
        set_data(d, 'phase_01/pore.test2', arr)
        set_data(d, 'pore.test3', arr)
        set_data(d, 'phase_02/pore.test3', arr)
        assert 'phase_02/pore.test3' in d.keys()

    def test_set_data_delete_array_with_None(self):
        d = {'pore.all': np.ones(10, dtype=bool)}
        set_data(d, 'pore.all', None)
        assert 'pore.all' not in d.keys()

    def test_set_data_delete_multiple_arrays_with_None(self):
        d = {'pore.all': np.ones(10, dtype=bool),
             'pore.test': np.ones(10)}
        set_data(d, '', {'pore.all': None, 'pore.test': None})
        assert 'pore.all' not in d.keys()
        assert 'pore.test' not in d.keys()

    def test_set_data_delete_group_with_None(self):
        d = {'network/pore.all': np.ones(10, dtype=bool)}
        set_data(d, 'network*', None)
        assert len(d.keys()) == 0

    def test_set_data_delete_some_arrays_with_wildcard_and_None(self):
        d = {'network/pore.foo': np.ones(10, dtype=bool),
             'network/pore.bar': np.ones(10, dtype=bool)}
        set_data(d, 'network/pore*', None)
        assert len(d.keys()) == 0

    def test_set_data_with_domain_notation(self):
        d = {'network/pore.all': np.ones(10, dtype=bool)}
        mask = np.zeros(10, dtype=bool)
        mask[[0, 1, 2, 3]] = True
        arr = np.zeros(10)
        d.update({'network/pore.left': mask})
        set_data(d, 'network/pore.test', arr)
        arr = np.ones(4)
        set_data(d, 'network/pore.test@left', arr)
        assert 'network/pore.test' in d.keys()
        assert np.sum(d['network/pore.test'] == 1) == 4

    def test_set_data_create_new_array_with_domain_notation(self):
        d = {'network/pore.all': np.ones(10, dtype=bool)}
        mask = np.zeros(10, dtype=bool)
        mask[[0, 1, 2, 3]] = True
        arr = np.zeros(4)
        d.update({'network/pore.left': mask})
        set_data(d, 'network/pore.test@left', arr)
        assert 'network/pore.test' in d.keys()
        assert np.sum(d['network/pore.test'] == 0) == 4

    def test_set_data_scalars_are_broadcast_to_full_arrays_all_present(self):
        d = {'network/pore.all': np.ones(10, dtype=bool)}
        set_data(d, 'network/pore.test', 1)
        assert d['network/pore.test'].shape == (10, )

    def test_set_data_scalars_are_blocked_if_size_unknown(self):
        d = {'network/pore.all': np.ones(10, dtype=bool)}
        # If value is a scalar throw error since actual size can't be inferred
        with pytest.raises(Exception):
            set_data(d, 'network/throat.test', 1)

    def test_set_data_vector_becomes_default_size(self):
        d = {'network/pore.all': np.ones(10, dtype=bool)}
        # If value is array-like, just write it AND create the 'all' array
        set_data(d, 'network/throat.test', [1, 1])
        assert d['network/throat.test'].shape[0] == 2
        with pytest.raises(Exception):
            set_data(d, 'network/throat.test', [1, 1, 1])

    def test_set_data_creating_deeply_nested_groups_on_fly(self):
        d = {'network/pore.all': np.ones(10, dtype=bool)}
        d2 = {'pore.all': np.ones(10, dtype=bool),
              'pore.test': 1.0,  # Try scalar!
              'param.MW': 29.1}  # Try param!
        set_data(d, 'network/phase1/phase2/phase3', d2)
        assert 'network/phase1/phase2/phase3/param.MW' in d.keys()

    def test_set_data_wrong_length_array_is_blocked(self):
        d = {'network/pore.all': np.ones(10, dtype=bool)}
        with pytest.raises(Exception):
            set_data(d, 'network/pore.test', np.ones(5))

    def test_set_data_valid_param(self):
        d = {'network/pore.all': np.ones(10, dtype=bool)}
        set_data(d, 'network/param.MW', 29.1)
        assert np.shape(d['network/param.MW']) == ()
        set_data(d, 'network/param.T', [273.1])
        assert np.shape(d['network/param.T']) == ()

    def test_set_data_invalid_param(self):
        d = {'network/pore.all': np.ones(10, dtype=bool)}
        with pytest.raises(Exception):
            set_data(d, 'network/param.MWs', [28.0, 32.0])

    def test_set_data_invalid_prefix(self):
        d = {'network/pore.all': np.ones(10, dtype=bool)}
        with pytest.raises(Exception):
            set_data(d, 'network/parameter.MWs', 29.1)

    def test_set_data_passing_dictionary_as_value_and_group_as_key(self):
        d = {'pore.all': np.ones(10, dtype=bool)}
        set_data(d, 'param.MW', 29.1)
        assert d['param.MW'] == 29.1
        set_data(d, '', {'attr.ID': 123, 'attr.parent': 321})
        assert 'attr.ID' in d.keys()
        assert 'attr.parent' in d.keys()

    def test_set_data_specifying_locs_on_existing_array(self):
        d = {'pore.all': np.ones(10, dtype=bool),
             'pore.test': np.ones(10, dtype=int)}
        set_data(d, 'pore.test', value=2, locs=[1, 2, 3])
        assert 'pore.test' in d.keys()
        assert (d['pore.test'] == 1).sum() == 7
        assert (d['pore.test'] == 2).sum() == 3
        assert np.all(d['pore.test'][[1, 2, 3]] == (2, 2, 2))

    def test_set_data_write_scalar_using_locs_on_nonexistant_array(self):
        d = {'pore.all': np.ones(10, dtype=bool)}
        set_data(d, 'pore.test', value=2, locs=[1, 2, 3])
        assert 'pore.test' in d.keys()
        assert np.isnan(d['pore.test']).sum() == 7
        assert (d['pore.test'] == 2.0).sum() == 3

    def test_set_data_write_array_using_locs_on_nonexistant_array(self):
        d = {'pore.all': np.ones(10, dtype=bool)}
        set_data(d, 'pore.test', value=[2, 2, 2], locs=[1, 2, 3])
        assert 'pore.test' in d.keys()
        assert np.isnan(d['pore.test']).sum() == 7
        assert (d['pore.test'] == 2.0).sum() == 3

    def test_set_data_write_conduit_prefix(self):
        d = {'throat.all': np.ones(10, dtype=bool)}
        arr = np.ones((10, 3))
        set_data(d, 'conduit.test', value=arr)

    def test_set_data_scalar_to_conduit(self):
        d = {'throat.all': np.ones(10, dtype=bool)}
        set_data(d, 'conduit.test', value=1)
        assert get_data(d, 'conduit.test').shape == (10, 3)

    def test_set_data_conduit_scalar_to_locs(self):
        d = {'throat.all': np.ones(10, dtype=bool)}
        set_data(d, 'conduit.test', value=1, locs=((1, 0), (2, 1)))
        assert get_data(d, 'conduit.test').shape == (10, 3)

    def test_set_data_param(self):
        d = {'pore.all': np.ones(10, dtype=bool)}
        set_data(d, 'param.test', value=1.0)
        assert 'param.test' in d.keys()
        assert np.size(d['param.test']) == 1
        with pytest.raises(Exception):
            set_data(d, 'param.test', value=[1.0, 2.0])

    def test_ensure_nested_dicts_get_flattened(self):
        d = {}
        e = {'phase': {'componentA': {'pore.temp': np.ones(10),
                                      'pore.blah': np.zeros(10)},
                       'componentB': {'pore.temp': np.ones(10),
                                      'pore.blah': np.zeros(10)}}}
        set_data(d, 'phase1', e)
        assert 'phase1/phase/componentA/pore.temp' in d.keys()
        assert 'phase1/phase/componentB/pore.temp' in d.keys()


if __name__ == '__main__':

    t = SetDataTests()
    t.run_all()
