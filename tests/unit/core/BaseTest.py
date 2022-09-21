import pytest
import numpy as np
import openpnm as op


class BaseTest:

    def setup_class(self):
        self.net = op.network.Cubic([3, 3, 3])
        self.net.set_label(throats=self.net.Ts, label='internal')
        self.net.set_label(pores=self.net.Ps, label='internal')

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_clear_model_data(self):
        pn = op.network.Cubic([5, 5, 5])
        phase = op.phase.Water(network=pn)
        a = len(phase)
        phase.clear(mode='models')
        assert len(phase) == (a - len(phase.models))
        # Clear non existing data
        phase.clear(mode='models')
        assert len(phase) == (a - len(phase.models))

    def test_clear_model_data_when_model_returns_dictionary(self):
        pn = op.network.Cubic([5, 5, 5])
        # geo = op.geometry.SpheresAndCylinders(network=pn, pores=pn.Ps, throats=pn.Ts)
        # geo.clear(mode='models')

    def test_pores(self):
        a = self.net.pores()
        assert np.all(a == np.arange(0, self.net.Np))

    def test_pores_one_label(self):
        a = self.net.pores(labels='top')
        assert np.all(a == [2, 5, 8, 11, 14, 17, 20, 23, 26])

    def test_pores_two_labels_or(self):
        a = self.net.pores(labels=['top', 'left'], mode='or')
        assert np.all(a == [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 14, 17, 20, 23, 26])

    def test_pores_two_labels_xnor(self):
        a = self.net.pores(labels=['top', 'left'], mode='xnor')
        assert np.all(a == [2, 5, 8])

    def test_pores_two_labels_not_xor(self):
        a = self.net.pores(labels=['top', 'left'], mode='xor')
        assert np.all(a == [0, 1, 3, 4, 6, 7, 11, 14, 17, 20, 23, 26])

    def test_pores_two_labels_nor(self):
        a = self.net.pores(labels=['top', 'left'], mode='nor')
        assert np.all(a == [9, 10, 12, 13, 15, 16, 18, 19, 21, 22, 24, 25])
        b = self.net.pores(labels=['top', 'left'], mode='or')
        c = self.net.to_mask(pores=a)*self.net.to_mask(pores=b)
        assert c.sum() == 0

    def test_pores_two_labels_nand(self):
        a = self.net.pores(labels=['top', 'left'], mode='nand')
        assert np.all(a == [0, 1, 3, 4, 6, 7, 11, 14, 17, 20, 23, 26])

    def test_pores_bad_mode(self):
        with pytest.raises(Exception):
            self.net.pores(labels=['top', 'front'], mode='bob')

    def test_pores_empty_list(self):
        a = self.net.pores(labels=[], mode='or')
        assert a.size == 0

    def test_pores_asmask(self):
        a = self.net.pores(labels=['top', 'front'], mode='or', asmask=True)
        assert a.sum() == 15
        b = self.net.pores(labels=['top', 'front'], mode='or')
        assert np.all(np.where(a)[0] == b)

    # def test_pores_and_throats_with_to_global(self):
    #     net = op.network.Cubic(shape=[2, 2, 2])
    #     geo1 = op.geometry.GenericGeometry(network=net,
    #                                        pores=[1, 3, 5, 7],
    #                                        throats=range(6))
    #     geo2 = op.geometry.GenericGeometry(network=net,
    #                                        pores=[0, 2, 4, 6],
    #                                        throats=range(6, 12))
    #     assert np.all(net.throats('geo_01') == geo1.throats(to_global=True))
    #     assert np.all(net.throats('geo_02') == geo2.throats(to_global=True))
    #     assert np.all(net.pores('geo_01') == geo1.pores(to_global=True))
    #     assert np.all(net.pores('geo_02') == geo2.pores(to_global=True))

    def test_throats(self):
        a = self.net.throats()
        assert np.all(a == np.arange(0, self.net.Nt))

    def test_throats_asmask(self):
        a = self.net.throats(labels=['all'], mode='or', asmask=True)
        assert a.sum() == 54
        b = self.net.throats(labels=['all'], mode='or')
        assert np.all(np.where(a)[0] == b)

    # def test_throats_one_label(self):
    #     self.geo['throat.label1'] = False
    #     self.geo['throat.label1'][[0, 1, 2, 3, 4, 5]] = True
    #     a = self.net.throats(labels='label1')
    #     assert np.all(a == [0, 1, 2, 3, 4, 5])
    #     self.geo.pop('throat.label1')

    # def test_throats_two_labels_or(self):
    #     self.geo['throat.label1'] = False
    #     self.geo['throat.label2'] = False
    #     self.geo['throat.label1'][[0, 1, 2, 3, 4, 5]] = True
    #     self.geo['throat.label2'][[6, 7, 8]] = True
    #     a = self.net.throats(labels=['label1', 'label2'], mode='or')
    #     assert np.all(a == [0, 1, 2, 3, 4, 5, 6, 7, 8])
    #     self.geo.pop('throat.label1')
    #     self.geo.pop('throat.label2')

    # def test_throats_two_labels_xnor(self):
    #     self.geo['throat.label1'] = False
    #     self.geo['throat.label2'] = False
    #     self.geo['throat.label1'][[0, 1, 2, 3, 4, 5]] = True
    #     self.geo['throat.label2'][[3, 4, 5, 6, 7, 8]] = True
    #     a = self.net.throats(labels=['label1', 'label2'], mode='xnor')
    #     assert np.all(a == [3, 4, 5])
    #     self.geo.pop('throat.label1')
    #     self.geo.pop('throat.label2')

    # def test_throats_two_labels_xor(self):
    #     self.geo['throat.label1'] = False
    #     self.geo['throat.label2'] = False
    #     self.geo['throat.label1'][[0, 1, 2, 3, 4, 5]] = True
    #     self.geo['throat.label2'][[3, 4, 5, 6, 7, 8]] = True
    #     a = self.net.throats(labels=['label1', 'label2'], mode='xor')
    #     assert np.all(a == [0, 1, 2, 6, 7, 8])
    #     self.geo.pop('throat.label1')
    #     self.geo.pop('throat.label2')

    def test_filter_by_label_pores_no_label(self):
        Ps = self.net.pores(['top', 'bottom', 'front'])
        with pytest.raises(Exception):
            self.net.filter_by_label(pores=Ps)

    def test_filter_by_label_pores_one_label_as_string(self):
        Ps = self.net.pores(['top', 'bottom', 'front'])
        a = self.net.filter_by_label(pores=Ps, labels='top', mode='or')
        b = [2, 5, 8, 11, 14, 17, 20, 23, 26]
        assert np.all(a == b)

    def test_filter_by_label_pores_one_label_as_list(self):
        Ps = self.net.pores(['top', 'bottom', 'front'])
        a = self.net.filter_by_label(pores=Ps, labels=['top'])
        b = [2, 5, 8, 11, 14, 17, 20, 23, 26]
        assert np.all(a == b)

    def test_filter_by_label_pores_two_labels_or(self):
        Ps = self.net.pores(['top', 'bottom', 'front'])
        a = self.net.filter_by_label(pores=Ps, labels=['top', 'bottom'],
                                     mode='or')
        b = [0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23, 24, 26]
        assert np.all(a == b)

    def test_filter_by_label_pores_two_labels_xnor(self):
        Ps = self.net.pores(['top', 'bottom', 'left'])
        a = self.net.filter_by_label(pores=Ps, labels=['top', 'left'],
                                     mode='xnor')
        b = [2, 5, 8]
        assert np.all(a == b)

    def test_filter_by_label_pores_two_labels_xnor_empty(self):
        Ps = self.net.pores(['top', 'bottom', 'left'])
        a = self.net.filter_by_label(pores=Ps, labels=['top', 'bottom'],
                                     mode='xnor')
        b = []
        assert np.all(a == b)

    def test_filter_by_label_pores_two_labels_xor(self):
        Ps = self.net.pores(['top', 'bottom', 'left'])
        a = self.net.filter_by_label(pores=Ps, labels=['top', 'left'],
                                     mode='xor')
        b = [0, 1, 3, 4, 6, 7, 11, 14, 17, 20, 23, 26]
        assert np.all(a == b)

    def test_filter_by_label_pores_two_labels_nor(self):
        Ps = self.net.pores(['top', 'bottom', 'left'])
        a = self.net.filter_by_label(pores=Ps, labels=['top', 'left'],
                                     mode='nor')
        b = [9, 12, 15, 18, 21, 24]
        assert np.all(a == b)

    def test_filter_by_label_empty_locations(self):
        a = self.net.filter_by_label(pores=[], labels='top')
        assert np.size(a) == 0

    def test_filter_by_label_pores_and_throats(self):
        with pytest.raises(Exception):
            self.net.filter_by_label(pores=[0, 1, 2], throats=[0, 1, 2])

    def test_tomask_pores(self):
        a = self.net.to_mask(pores=self.net.pores('top'))
        assert np.sum(a) == 9

    # def test_tomask_throats(self):
    #     self.geo['throat.label1'] = False
    #     self.geo['throat.label1'][[0, 1, 2, 3, 4, 5]] = True
    #     a = self.net.to_mask(throats=self.net.throats('label1'))
    #     assert np.sum(a) == 6
    #     self.geo.pop('throat.label1')

    # def test_tomask_pores_and_throats(self):
    #     with pytest.raises(Exception):
    #         _ = self.net.to_mask(throats=[0, 1, 2], pores=[0, 1, 2])

    def test_toindices_pores(self):
        mask = np.zeros((self.net.Np), dtype=bool)
        Ps = [0, 3, 6]
        mask[Ps] = True
        a = self.net.to_indices(mask)
        assert np.all(a == Ps)

    def test_toindices_throats(self):
        mask = np.zeros((self.net.Nt), dtype=bool)
        Ts = [0, 3, 6]
        mask[Ts] = True
        a = self.net.to_indices(mask)
        assert np.all(a == Ts)

    def test_toindices_float_mask(self):
        mask = (np.random.rand(self.net.Np) < 0.5)
        inds_in = np.where(mask)[0]
        inds_out = self.net.to_indices(mask*1.0)
        assert np.all(inds_in == inds_out)

    # def test_toindices_invalid_mask(self):
    #     mask = self.net.Np
    #     with pytest.raises(Exception):
    #         self.net.to_indices(mask)

    # def test_toindices_wrong_mask(self):
    #     mask = np.zeros((self.net.Nt)-2, dtype=bool)
    #     mask[[0, 3, 6]] = True
    #     with pytest.raises(Exception):
    #         self.net.to_indices(mask)

    def test_count(self):
        with pytest.raises(Exception):
            self.net._count()

    def test_num_pores(self):
        a = self.net.num_pores()
        assert a == 27

    def test_num_pores_one_label(self):
        a = self.net.num_pores(labels='top')
        assert a == 9

    def test_num_pores_two_labels_or(self):
        a = self.net.num_pores(labels=['top', 'front'], mode='or')
        assert a == 15

    def test_num_pores_two_labels_xnor(self):
        a = self.net.num_pores(labels=['top', 'front'], mode='xnor')
        assert a == 3

    def test_num_pores_two_labels_xor(self):
        a = self.net.num_pores(labels=['top', 'front'],
                               mode='xor')
        assert a == 12

    def test_num_pores_two_labels_nor(self):
        a = self.net.num_pores(labels=['top', 'front'], mode='nor')
        assert a == 12

    def test_num_throats(self):
        a = self.net.num_throats()
        assert a == 54

    # def test_num_throats_one_label(self):
    #     self.geo['throat.label1'] = False
    #     self.geo['throat.label1'][[0, 1, 2, 3, 4, 5]] = True
    #     a = self.net.num_throats(labels='label1')
    #     assert a == 6
    #     self.geo.pop('throat.label1')

    # def test_num_throats_two_labels_or(self):
    #     self.geo['throat.label1'] = False
    #     self.geo['throat.label2'] = False
    #     self.geo['throat.label1'][[0, 1, 2, 3, 4, 5]] = True
    #     self.geo['throat.label2'][[3, 4, 5, 6, 7, 8]] = True
    #     a = self.net.num_throats(labels=['label1', 'label2'], mode='or')
    #     assert a == 9
    #     self.geo.pop('throat.label1')
    #     self.geo.pop('throat.label2')

    # def test_num_throats_two_labels_xnor(self):
    #     self.geo['throat.label1'] = False
    #     self.geo['throat.label2'] = False
    #     self.geo['throat.label1'][[0, 1, 2, 3, 4, 5]] = True
    #     self.geo['throat.label2'][[3, 4, 5, 6, 7, 8]] = True
    #     a = self.net.num_throats(labels=['label1', 'label2'],
    #                              mode='xnor')
    #     assert a == 3
    #     self.geo.pop('throat.label1')
    #     self.geo.pop('throat.label2')

    # def test_num_throats_two_labels_xor(self):
    #     self.geo['throat.label1'] = False
    #     self.geo['throat.label2'] = False
    #     self.geo['throat.label1'][[0, 1, 2, 3, 4, 5]] = True
    #     self.geo['throat.label2'][[3, 4, 5, 6, 7, 8]] = True
    #     a = self.net.num_throats(labels=['label1', 'label2'],
    #                              mode='xor')
    #     assert a == 6
    #     self.geo.pop('throat.label1')
    #     self.geo.pop('throat.label2')

    # def test_num_throats_two_labels_nor(self):
    #     self.geo['throat.label1'] = False
    #     self.geo['throat.label2'] = False
    #     self.geo['throat.label1'][[0, 1, 2, 3, 4, 5]] = True
    #     self.geo['throat.label2'][[3, 4, 5, 6, 7, 8]] = True
    #     a = self.net.num_throats(labels=['label1', 'label2'],
    #                              mode='nor')
    #     assert a == 45
    #     self.geo.pop('throat.label1')
    #     self.geo.pop('throat.label2')

    def test_keys_mode_skip(self):
        a = self.net.keys()
        assert 'dict_keys' in str(type(a))

    def test_keys_mode_props(self):
        a = self.net.keys(mode='props')
        assert 'dict_keys' not in str(type(a))
        b = [i for i in a if self.net[i].dtype != bool]
        assert a == b

    def test_keys_mode_labels(self):
        a = self.net.keys(mode='labels')
        assert 'dict_keys' not in str(type(a))
        b = [i for i in a if self.net[i].dtype == bool]
        assert a == b

    # def test_keys_element_pores_mode_all(self):
    #     a = self.net.keys(element='pores', mode='all')
    #     b = [i.split('.')[0] for i in a]
    #     assert set(b) == {'pore'}

    # def test_keys_element_throats_mode_all(self):
    #     a = self.net.keys(element='throats', mode='all')
    #     b = [i.split('.')[0] for i in a]
    #     assert set(b) == {'throat'}

    def test_keys_mode_props_and_labels(self):
        a = self.net.keys(mode=['props', 'labels'])
        b = list(self.net.keys())
        assert set(a) == set(b)

    # def test_props_all(self):
    #     a = self.geo.props()
    #     b = ['pore.diameter', 'pore.volume', 'throat.cross_sectional_area',
    #          'throat.diameter']
    #     assert sorted(a) == b

    # def test_props_models(self):
    #     a = self.geo.props(mode='models')
    #     b = ['pore.volume', 'throat.cross_sectional_area']
    #     assert sorted(a) == sorted(b)

    # def test_props_constants(self):
    #     a = self.geo.props(mode='constants')
    #     b = ['pore.diameter', 'throat.diameter']
    #     assert sorted(a) == sorted(b)

    # def test_props_pores_all(self):
    #     a = self.geo.props(element='pores')
    #     b = ['pore.diameter', 'pore.volume']
    #     assert sorted(a) == sorted(b)

    # def test_props_pores_models(self):
    #     a = self.geo.props(element='pores', mode='models')
    #     b = ['pore.volume']
    #     assert sorted(a) == sorted(b)

    # def test_props_pores_constants(self):
    #     a = self.geo.props(element='pores', mode='constants')
    #     b = ['pore.diameter']
    #     assert sorted(a) == sorted(b)

    def test_props_hidden_keys(self):
        self.net['pore._blah'] = 1.0
        assert 'pore._blah' not in self.net.__str__()
        assert 'pore._blah' not in self.net.props()
        assert 'pore._blah' in self.net.keys()

    def test_labels(self):
        a = self.net.labels()
        assert 'pore.top' in a

    def test_labels_on_pores(self):
        a = self.net.labels(element='pore')
        b = ['pore.back', 'pore.bottom', 'pore.front',
             'pore.internal', 'pore.left', 'pore.right',
             'pore.top', 'pore.surface']
        assert sorted(a) == sorted(b)

    def test_labels_on_throats(self):
        a = self.net.labels(element='throat')
        b = ['throat.internal', 'throat.surface']
        assert sorted(a) == sorted(b)

    def test_labels_on_pores_and_throats(self):
        with pytest.raises(Exception):
            self.net.labels(pores=[0, 1], throats=[0, 1])

    # def test_labels_on_foo(self):
    #     with pytest.raises(Exception):
    #         self.net.labels(element='foo')

    def test_labels_on_all_pores(self):
        a = self.net.labels(pores=self.net.Ps)
        b = ['pore.back', 'pore.bottom', 'pore.front',
             'pore.internal', 'pore.left', 'pore.right','pore.top',
             'pore.surface']
        assert sorted(a) == sorted(b)

    def test_labels_on_all_throats(self):
        a = self.net.labels(throats=self.net.Ts)
        b = ['throat.internal', 'throat.surface']
        assert sorted(a) == sorted(b)

    def test_labels_on_one_pore(self):
        a = self.net.labels(pores=0)
        b = ['pore.bottom', 'pore.front',
             'pore.internal', 'pore.surface',
             'pore.left']
        assert sorted(a) == sorted(b)

    def test_labels_on_list_of_pores(self):
        a = self.net.labels(pores=[0, 1])
        b = ['pore.bottom', 'pore.front',
             'pore.internal', 'pore.surface',
             'pore.left']
        assert sorted(a) == sorted(b)

    def test_labels_locations_boolean(self):
        ind = np.zeros((self.net.Np), dtype=bool)
        ind[[0, 1]] = True
        a = self.net.labels(pores=ind)
        b = ['pore.bottom', 'pore.front',
             'pore.internal', 'pore.surface',
             'pore.left']
        assert sorted(a) == sorted(b)

    def test_labels_pores_mode_or(self):
        a = self.net.labels(pores=[0, 1, 2], mode='or')
        b = ['pore.bottom', 'pore.front',
             'pore.internal', 'pore.surface',
             'pore.left', 'pore.top']
        assert sorted(a) == sorted(b)

    def test_labels_pores_mode_and(self):
        a = self.net.labels(pores=[0, 1, 2], mode='and')
        b = ['pore.front', 'pore.internal', 'pore.left',
             'pore.surface']
        assert sorted(a) == sorted(b)

    def test_labels_pores_mode_xor(self):
        a = self.net.labels(pores=[0, 1, 2], mode='xor')
        b = ['pore.bottom', 'pore.top']
        assert sorted(a) == sorted(b)

    def test_labels_pores_mode_nand(self):
        a = self.net.labels(pores=[0, 1, 2], mode='nand')
        b = []
        assert sorted(a) == sorted(b)

    def test_labels_pores_mode_xnor(self):
        a = self.net.labels(pores=[0, 1, 2], mode='xnor')
        b = ['pore.front', 'pore.internal',
             'pore.surface', 'pore.left']
        assert sorted(a) == sorted(b)

    def test_labels_pores_mode_nor(self):
        a = self.net.labels(pores=[0, 1, 2], mode='nor')
        b = ['pore.back', 'pore.right']
        assert sorted(a) == sorted(b)

    def test_labels_pores_mode_foo(self):
        with pytest.raises(Exception):
            self.net.labels(pores=[0, 1], mode='foo')

    def test_labels_hidden_key(self):
        self.net['pore._foo'] = True
        assert 'pore._foo' not in self.net.__str__()
        assert 'pore._foo' in self.net.keys()

    def test_parse_indices_boolean(self):
        b = np.array([True, True, True])
        with pytest.raises(Exception):
            self.net._parse_indices(indices=b)
        b = np.zeros((self.net.Np,), dtype=bool)
        assert len(self.net._parse_indices(indices=b)) == 0
        b = np.zeros((self.net.Nt,), dtype=bool)
        b[[0, 1, 2]] = True
        assert np.shape(self.net._parse_indices(indices=b)) == (3,)

    def test_parse_indices_None(self):
        assert len(self.net._parse_indices(indices=None)) == 0

    def test_parse_indices_int(self):
        a = self.net._parse_indices(indices=0)
        assert isinstance(a, np.ndarray)
        assert np.all(a == 0)

    def test_parse_indices_list(self):
        a = self.net._parse_indices(indices=[0, 1])
        assert isinstance(a, np.ndarray)
        assert np.all(a == [0, 1])

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
        a = self.net._parse_mode(mode='or')
        assert a == ['or']

    def test_parse_mode_single(self):
        a = self.net._parse_mode(mode=['or', 'xnor'])
        assert sorted(a) == ['or', 'xnor']
        with pytest.raises(Exception):
            a = self.net._parse_mode(mode=['or1', 'or2'], single=True)
        a = self.net._parse_mode(mode=['or1'], single=True)
        assert a == 'or1'

    def test_parse_mode_allowed(self):
        allowed = ['a', 'b', 'c']
        with pytest.raises(Exception):
            self.net._parse_mode(mode=['a', 'd'], allowed=allowed)

    def test_parse_mode_duplicate(self):
        a = self.net._parse_mode(mode=['or', 'or'])
        assert a == ['or']
        a = self.net._parse_mode(mode=['or', 'or'], single=True)
        assert a == 'or'

    def test_setitem_wrong_prefix(self):
        with pytest.raises(Exception):
            self.geo['pore2.test'] = 0

    def test_setitem_wrong_length(self):
        with pytest.raises(Exception):
            self.net['pore.test'] = np.ones((self.net.Np+1))
        assert 'pore.test' not in self.net.keys()

    # def test_setitem_replace_all(self):
    #     array_len = np.size(self.geo['pore.all'])
    #     self.geo['pore.all'] = np.ones((self.geo.Np+1))
    #     assert np.size(self.geo['pore.all']) == array_len

    def test_setitem_overwrite_into_all(self):
        pass
        # This test will fail as there is currently no way to prevent this
        # array_sum = np.sum(self.geo['pore.all'])
        # self.geo['pore.all'][0] = False
        # assert np.sum(self.geo['pore.all']) == array_sum

    # def test_setitem_subdict_conflicts(self):
    #     self.geo['pore.foo'] = 1
    #     with pytest.raises(Exception):
    #         self.geo['pore.foo.bar'] = 1
    #     self.geo['throat.foo.bar'] = 1
    #     with pytest.raises(Exception):
    #         self.geo['throat.foo'] = 1

    # def test_object_name_name_conflict(self):
    #     with pytest.raises(Exception):
    #         self.geo.name = self.net.name

    # def test_object_name_array_conflict(self):
    #     with pytest.raises(Exception):
    #         self.geo.name = 'coords'
    #     Np = self.geo.Np
    #     Nt = self.geo.Nt
    #     assert self.geo.Np == Np
    #     assert self.geo.Nt == Nt

    def test_get_indices_wildcard(self):
        a = self.net._get_indices(element='pore', labels='ba*')
        assert np.all(a == [6, 7, 8, 15, 16, 17, 24, 25, 26])
        b = self.net._get_indices(element='pore', labels='*ck')
        assert np.all(a == b)

    def test_write_dict(self):
        self.net['pore.test_dict'] = {'test1': 1, 'test2': self.net.Ps}
        assert 'pore.test_dict.test1' in self.net.keys()
        assert self.net['pore.test_dict.test1'].shape == (self.net.Np, )
        assert 'pore.test_dict.test2' in self.net.keys()

    # def test_map_pores(self):
    #     a = self.geo21['pore._id']
    #     b = self.geo22['pore._id']
    #     assert a.size == self.geo21.Np
    #     assert b.size == self.geo22.Np
    #     assert ~np.any(np.in1d(a, b))
    #     Pgeo21 = self.net2.map_pores(pores=self.geo21.Ps, origin=self.geo21)
    #     assert np.all(Pgeo21 == self.net2.pores(self.geo21.name))
    #     Pgeo22 = self.net2.map_pores(pores=self.geo22.Ps, origin=self.geo22)
    #     assert np.all(Pgeo22 == self.net2.pores(self.geo22.name))

    # def test_map_throats(self):
    #     a = self.geo21['throat._id']
    #     assert a.size == self.geo21.Nt
    #     Tgeo21 = self.net2.map_throats(throats=self.geo21.Ts,
    #                                    origin=self.geo21)
    #     assert np.all(Tgeo21 == self.net2.throats(self.geo21.name))

    # def test_map_pores_unfiltered(self):
    #     b = self.net.map_pores(pores=self.geo.Ps, origin=self.geo,
    #                            filtered=False)
    #     assert np.all(b.indices == self.net.pores(self.geo.name))
    #     assert b.mask.size == self.geo.Np

    # def test_map_pores_unfiltered_missing(self):
    #     Ps = self.net2.Ps[15:20]
    #     b = self.geo22.map_pores(pores=Ps, origin=self.net2, filtered=False)
    #     assert sum(b.mask) == 2
    #     assert len(b.mask) == 5

    # def test_map_pores_reverse(self):
    #     Ps = self.net2.Ps[:5]
    #     b = self.geo21.map_pores(pores=Ps, origin=self.net2)
    #     assert np.all(b == [0, 1, 2, 3, 4])
    #     Ps = self.net2.Ps[-5:]
    #     b = self.geo22.map_pores(pores=Ps, origin=self.net2)
    #     assert np.all(b == [4, 5, 6, 7, 8])

    # def test_map_pores_missing(self):
    #     Ps = self.net2.Ps[:5]
    #     b = self.geo22.map_pores(pores=Ps, origin=self.net2)
    #     assert len(b) == 0

    # def test_getitem_with_no_matches(self):
    #     self.geo.pop('pore.blah', None)
    #     with pytest.raises(KeyError):
    #         _ = self.geo['pore.blah']

    def test_getitem_with_numeric_value(self):
        assert self.net[44] == 44
        assert self.net[44.4] == 44.4
        assert self.net[4j] == 4j

    # def test_interpolate_data(self):
    #     self.geo['throat.tester'] = np.linspace(0, 1.0, self.geo.network.Nt)
    #     self.geo['pore.tester'] = np.linspace(0, 1.0, self.geo.network.Np)
    #     a = self.geo.interpolate_data(propname='throat.tester')
    #     assert a.size == self.geo.Np
    #     assert np.isclose(a.mean(), 0.5)
    #     a = self.geo.interpolate_data(propname='pore.tester')
    #     assert a.size == self.geo.Nt
    #     assert np.isclose(a.mean(), 0.5)

    def test_get_string(self):
        a = self.net.get('pore.coords')
        assert a.shape == (self.net.Np, 3)

    def test_subdict_getitem_on_network_from_network(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        pn['pore.foo.bar'] = 1
        pn['pore.foo.baz'] = 2
        d = pn['pore.foo']
        assert len(d) == 2
        assert 'bar' in d.keys()
        assert 'baz' in d.keys()

    def test_subdict_lookup_errors(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        pn['pore.foo.bar'] = 1
        pn['pore.foo.baz'] = 2
        with pytest.raises(KeyError):
            _ = pn['pore.foo.b']
        with pytest.raises(KeyError):
            _ = pn['pore.fo']

    def test_set_label_add_to_pores(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        pn.set_label(label='tester', pores=[1, 2])
        assert pn['pore.tester'].sum() == 2

    def test_set_label_add_to_throats(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        pn.set_label(label='tester', throats=[1, 2])
        assert pn['throat.tester'].sum() == 2

    def test_set_label_overwrite_on_pores(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        pn.set_label(label='tester', pores=[1, 2])
        pn.set_label(label='tester', pores=[2, 3, 4], mode='overwrite')
        assert pn['pore.tester'].sum() == 3

    def test_set_label_overwrite_on_throats(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        pn.set_label(label='tester', throats=[1, 2])
        pn.set_label(label='tester', throats=[2, 3, 4], mode='overwrite')
        assert pn['throat.tester'].sum() == 3

    def test_set_label_remove_from_pores(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        pn.set_label(label='tester', pores=[1, 2])
        assert pn['pore.tester'].sum() == 2
        pn.set_label(label='tester', pores=[1, 2, 3], mode='remove')
        assert pn['pore.tester'].sum() == 0

    def test_set_label_remove_from_throats(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        pn.set_label(label='tester', throats=[1, 2])
        assert pn['throat.tester'].sum() == 2
        pn.set_label(label='tester', throats=[1, 2, 3], mode='remove')
        assert pn['throat.tester'].sum() == 0

    def test_set_label_purge_from_pores(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        pn.set_label(label='tester', pores=[1, 2])
        assert pn['pore.tester'].sum() == 2
        pn.set_label(label='tester', mode='purge')
        assert 'pore.tester' not in pn.keys()

    def test_set_label_purge_from_throats(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        pn.set_label(label='tester', throats=[1, 2])
        assert pn['throat.tester'].sum() == 2
        pn.set_label(label='tester', mode='purge')
        assert 'throat.tester' not in pn.keys()

    def test_set_label_purge_nonexistent_label(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        pn.set_label(label='tester', mode='purge')
        # Should only issue warning

    def test_renaming_to_current_name_is_allowed(self):
        obj = op.core.Base2(name="temp")
        obj.name = "temp"

    def test_get_conduit_data(self):
        pn = op.network.Cubic(shape=[3, 3, 3])
        pn['pore.diameter'] = 2.0
        pn['throat.diameter'] = 2.0
        a = pn.get_conduit_data('diameter')
        assert np.all(a.shape == (54, 3))
        assert a.sum() == np.prod(a.shape)*2.0
        pn['pore.volume'] = 2.0
        pn['throat.volume'] = 2.0
        b = pn.get_conduit_data('volume')
        assert np.all(b.shape == (54, 3))
        assert b.sum() == 2*np.prod(b.shape)
        with pytest.raises(KeyError):
            pn.get_conduit_data('blah')


if __name__ == '__main__':

    t = BaseTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
