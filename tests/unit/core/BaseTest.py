import openpnm as op
import scipy as sp
import pytest


class BaseTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.diameter'] = sp.rand(self.net.Np)
        self.geo.add_model(propname='pore.volume',
                           model=op.models.geometry.pore_volume.sphere)
        self.geo['throat.diameter'] = sp.rand(self.net.Nt)
        self.geo.add_model(propname='throat.area',
                           model=op.models.geometry.throat_area.cylinder)
        self.geo.regenerate_models()
        self.geo['throat.label1'] = False
        self.geo['throat.label2'] = False
        self.geo['throat.label1'][0:6] = True
        self.geo['throat.label2'][3:9] = True
        self.net1 = op.network.Cubic(shape=[3, 3, 3])
        self.geo1 = op.geometry.GenericGeometry(network=self.net1,
                                                pores=self.net1.Ps,
                                                throats=self.net1.Ts)
        self.phase1 = op.phases.GenericPhase(network=self.net1)
        self.phase2 = op.phases.GenericPhase(network=self.net1)
        self.phys1 = op.physics.GenericPhysics(network=self.net1,
                                               geometry=self.geo1,
                                               phase=self.phase1)
        self.phys2 = op.physics.GenericPhysics(network=self.net1,
                                               geometry=self.geo1,
                                               phase=self.phase2)
        self.net2 = op.network.Cubic(shape=[3, 3, 3])
        Ps = sp.arange(0, 18)
        Ts = self.net2.find_neighbor_pores(Ps, mode='union')
        self.geo21 = op.geometry.GenericGeometry(network=self.net2,
                                                 pores=Ps,
                                                 throats=Ts)
        Ps = sp.arange(18, 27)
        Ts = self.net2.find_neighbor_pores(Ps, mode='intersection')
        self.geo22 = op.geometry.GenericGeometry(network=self.net2,
                                                 pores=Ps,
                                                 throats=Ts)

    def teardown_class(self):
        ws = op.core.Workspace()
        ws.clear()

    def test_pores(self):
        a = self.net.pores()
        assert sp.all(a == sp.arange(0, self.net.Np))

    def test_pores_one_label(self):
        a = self.net.pores(labels='top')
        assert sp.all(a == [2, 5, 8, 11, 14, 17, 20, 23, 26])

    def test_pores_two_labels_union(self):
        a = self.net.pores(labels=['top', 'front'], mode='union')
        assert sp.all(a == [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 14, 17, 20, 23, 26])

    def test_pores_two_labels_intersection(self):
        a = self.net.pores(labels=['top', 'front'], mode='intersection')
        assert sp.all(a == [2, 5, 8])

#    def test_pores_two_labels_not_intersection(self):
#        a = self.net.pores(labels=['top', 'front'], mode='not_intersection')
#        assert sp.all(a == [0, 1, 3, 4, 6, 7, 11, 14, 17, 20, 23, 26])

    def test_pores_two_labels_difference(self):
        a = self.net.pores(labels=['top', 'front'], mode='difference')
        assert sp.all(a == [9, 10, 12, 13, 15, 16, 18, 19, 21, 22, 24, 25])

    def test_throats(self):
        a = self.net.throats()
        assert sp.all(a == sp.arange(0, self.net.Nt))

    def test_throats_one_label(self):
        a = self.net.throats(labels='label1')
        assert sp.all(a == [0, 1, 2, 3, 4, 5])

    def test_throats_two_labels_union(self):
        a = self.net.throats(labels=['label1', 'label2'], mode='union')
        assert sp.all(a == [0, 1, 2, 3, 4, 5, 6, 7, 8])

    def test_throats_two_labels_intersection(self):
        a = self.net.throats(labels=['label1', 'label2'], mode='intersection')
        assert sp.all(a == [3, 4, 5])

#    def test_throats_two_labels_not_intersection(self):
#        a = self.net.throats(labels=['label1', 'label2'],
#                             mode='not_intersection')
#        assert sp.all(a == [0, 1, 2, 6, 7, 8])

    def test_filter_by_label_pores_no_label(self):
        Ps = self.net.pores(['top', 'bottom', 'front'])
        with pytest.raises(Exception):
            self.net.filter_by_label(pores=Ps)

    def test_filter_by_label_pores_one_label_as_string(self):
        Ps = self.net.pores(['top', 'bottom', 'front'])
        a = self.net.filter_by_label(pores=Ps, labels='top')
        b = [2, 5, 8, 11, 14, 17, 20, 23, 26]
        assert sp.all(a == b)

    def test_filter_by_label_pores_one_label_as_list(self):
        Ps = self.net.pores(['top', 'bottom', 'front'])
        a = self.net.filter_by_label(pores=Ps, labels=['top'])
        b = [2, 5, 8, 11, 14, 17, 20, 23, 26]
        assert sp.all(a == b)

    def test_filter_by_label_pores_two_labels_union(self):
        Ps = self.net.pores(['top', 'bottom', 'front'])
        a = self.net.filter_by_label(pores=Ps,
                                     labels=['top', 'bottom'],
                                     mode='union')
        b = [0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23, 24, 26]
        assert sp.all(a == b)

    def test_filter_by_label_pores_two_labels_intersection(self):
        Ps = self.net.pores(['top', 'bottom', 'front'])
        a = self.net.filter_by_label(pores=Ps,
                                     labels=['top', 'front'],
                                     mode='intersection')
        b = [2, 5, 8]
        assert sp.all(a == b)

    def test_filter_by_label_pores_two_labels_intersection_empty(self):
        Ps = self.net.pores(['top', 'bottom', 'front'])
        a = self.net.filter_by_label(pores=Ps,
                                     labels=['top', 'bottom'],
                                     mode='intersection')
        b = []
        assert sp.all(a == b)

#    def test_filter_by_label_pores_two_labels_not_intersection(self):
#        Ps = self.net.pores(['top', 'bottom', 'front'])
#        a = self.net.filter_by_label(pores=Ps,
#                                     labels=['top', 'front'],
#                                     mode='not_intersection')
#        b = [0, 1, 3, 4, 6, 7, 11, 14, 17, 20, 23, 26]
#        assert sp.all(a == b)

    def test_filter_by_label_pores_two_labels_complement(self):
        Ps = self.net.pores(['top', 'bottom', 'front'])
        a = self.net.filter_by_label(pores=Ps,
                                     labels=['top', 'front'],
                                     mode='complement')
        b = [9, 12, 15, 18, 21, 24]
        assert sp.all(a == b)

    def test_filter_by_label_empty_locations(self):
        a = self.net.filter_by_label(pores=[], labels='top')
        assert sp.size(a) == 0

    def test_tomask_pores(self):
        a = self.net.tomask(pores=self.net.pores('top'))
        assert sp.sum(a) == 9

    def test_tomask_throats(self):
        a = self.net.tomask(throats=self.net.throats('label1'))
        assert sp.sum(a) == 6

    def test_toindices_pores(self):
        mask = sp.zeros((self.net.Np), dtype=bool)
        Ps = [0, 3, 6]
        mask[Ps] = True
        a = self.net.toindices(mask)
        assert sp.all(a == Ps)

    def test_toindices_throats(self):
        mask = sp.zeros((self.net.Nt), dtype=bool)
        Ts = [0, 3, 6]
        mask[Ts] = True
        a = self.net.toindices(mask)
        assert sp.all(a == Ts)

    def test_toindices_float_mask(self):
        mask = (sp.rand(self.net.Np) < 0.5)
        inds_in = sp.where(mask)[0]
        inds_out = self.net.toindices(mask*1.0)
        assert sp.all(inds_in == inds_out)

    def test_toindices_invalid_mask(self):
        mask = self.net.Np
        with pytest.raises(Exception):
            self.net.toindices(mask)

    def test_toindices_wrong_mask(self):
        mask = sp.zeros((self.net.Nt)-2, dtype=bool)
        mask[[0, 3, 6]] = True
        with pytest.raises(Exception):
            self.net.toindices(mask)

    def test_count(self):
        with pytest.raises(Exception):
            self.net._count()

    def test_num_pores(self):
        a = self.net.num_pores()
        assert a == 27

    def test_num_pores_one_label(self):
        a = self.net.num_pores(labels='top')
        assert a == 9

    def test_num_pores_two_labels_union(self):
        a = self.net.num_pores(labels=['top', 'front'], mode='union')
        assert a == 15

    def test_num_pores_two_labels_intersection(self):
        a = self.net.num_pores(labels=['top', 'front'], mode='intersection')
        assert a == 3

#    def test_num_pores_two_labels_not_intersection(self):
#        a = self.net.num_pores(labels=['top', 'front'],
#                               mode='not_intersection')
#        assert a == 12

    def test_num_pores_two_labels_difference(self):
        a = self.net.num_pores(labels=['top', 'front'], mode='difference')
        assert a == 12

    def test_num_throats(self):
        a = self.net.num_throats()
        assert a == 54

    def test_num_throats_one_label(self):
        a = self.net.num_throats(labels='label1')
        assert a == 6

    def test_num_throats_two_labels_union(self):
        a = self.net.num_throats(labels=['label1', 'label2'], mode='union')
        assert a == 9

    def test_num_throats_two_labels_intersection(self):
        a = self.net.num_throats(labels=['label1', 'label2'],
                                 mode='intersection')
        assert a == 3

#    def test_num_throats_two_labels_notintersection(self):
#        a = self.net.num_throats(labels=['label1', 'label2'],
#                                 mode='not_intersection')
#        assert a == 6

    def test_num_throats_two_labels_difference(self):
        a = self.net.num_throats(labels=['label1', 'label2'],
                                 mode='difference')
        assert a == 45

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

    def test_keys_element_pores_mode_all(self):
        a = self.net.keys(element='pores', mode='all')
        b = [i.split('.')[0] for i in a]
        assert set(b) == {'pore'}

    def test_keys_element_throats_mode_all(self):
        a = self.net.keys(element='throats', mode='all')
        b = [i.split('.')[0] for i in a]
        assert set(b) == {'throat'}

    def test_keys_mode_props_and_labels(self):
        a = self.net.keys(mode=['props', 'labels'])
        b = list(self.net.keys())
        assert set(a) == set(b)

    def test_props_all(self):
        a = self.geo.props()
        assert sorted(a) == ['pore.diameter', 'pore.volume',
                             'throat.area', 'throat.diameter']

    def test_props_models(self):
        a = self.geo.props(mode='models')
        b = ['pore.volume', 'throat.area']
        assert sorted(a) == sorted(b)

    def test_props_constants(self):
        a = self.geo.props(mode='constants')
        b = ['pore.diameter', 'throat.diameter']
        assert sorted(a) == sorted(b)

    def test_props_pores_all(self):
        a = self.geo.props(element='pores')
        b = ['pore.diameter', 'pore.volume']
        assert sorted(a) == sorted(b)

    def test_props_pores_models(self):
        a = self.geo.props(element='pores', mode='models')
        b = ['pore.volume']
        assert sorted(a) == sorted(b)

    def test_props_pores_constants(self):
        a = self.geo.props(element='pores', mode='constants')
        b = ['pore.diameter']
        assert sorted(a) == sorted(b)

    def test_props_hidden_keys(self):
        self.net['pore._blah'] = 1.0
        assert 'pore._blah' not in self.net.__str__()
        assert 'pore._blah' not in self.net.props()
        assert 'pore._blah' in self.net.keys()

    def test_labels(self):
        a = self.net.labels()
        assert 'pore.top' in a

    def test_labels_on_pores(self):
        a = self.net.labels(element='pores')
        b = ['pore.all', 'pore.back', 'pore.bottom', 'pore.front',
             'pore.internal', 'pore.left', 'pore.right', 'pore.'+self.geo.name,
             'pore.top']
        assert sorted(a) == sorted(b)

    def test_labels_on_throats(self):
        a = self.net.labels(element='throats')
        b = ['throat.all', 'throat.internal', 'throat.'+self.geo.name]
        assert sorted(a) == sorted(b)

    def test_labels_on_foo(self):
        with pytest.raises(Exception):
            self.net.labels(element='foo')

    def test_labels_on_all_pores(self):
        a = self.net.labels(pores=self.net.Ps)
        b = ['pore.all', 'pore.back', 'pore.bottom', 'pore.front',
             'pore.internal', 'pore.left', 'pore.right', 'pore.'+self.geo.name,
             'pore.top']
        assert sorted(a) == sorted(b)

    def test_labels_on_all_throats(self):
        a = self.net.labels(throats=self.net.Ts)
        b = ['throat.all', 'throat.internal', 'throat.'+self.geo.name]
        assert sorted(a) == sorted(b)

    def test_labels_on_one_pore(self):
        a = self.net.labels(pores=0)
        b = ['pore.all', 'pore.bottom', 'pore.front', 'pore.internal',
             'pore.left', 'pore.'+self.geo.name]
        assert sorted(a) == sorted(b)

    def test_labels_on_list_of_pores(self):
        a = self.net.labels(pores=[0, 1])
        b = ['pore.all', 'pore.bottom', 'pore.front', 'pore.internal',
             'pore.left', 'pore.'+self.geo.name]
        assert sorted(a) == sorted(b)

    def test_labels_locations_boolean(self):
        ind = sp.zeros((self.net.Np), dtype=bool)
        ind[[0, 1]] = True
        a = self.net.labels(pores=ind)
        b = ['pore.all', 'pore.bottom', 'pore.front', 'pore.internal',
             'pore.left', 'pore.'+self.geo.name]
        assert sorted(a) == sorted(b)

    def test_labels_pores_mode_union(self):
        a = self.net.labels(pores=[0, 1, 2], mode='union')
        b = ['pore.all', 'pore.bottom', 'pore.front', 'pore.internal',
             'pore.left', 'pore.'+self.geo.name, 'pore.top']
        assert sorted(a) == sorted(b)

    def test_labels_pores_mode_intersection(self):
        a = self.net.labels(pores=[0, 1, 2], mode='intersection')
        b = ['pore.all', 'pore.front', 'pore.internal', 'pore.left',
             'pore.'+self.geo.name]
        assert sorted(a) == sorted(b)

    def test_labels_pores_mode_count(self):
        a = self.net.labels(pores=[0, 1, 2], mode='count')
        assert sp.all(a == [6, 5, 6])

    def test_labels_pores_mode_mask(self):
        a = self.net.labels(pores=[0, 1], mode='mask')
        assert sp.sum(a) == 11

    def test_labels_pores_mode_difference(self):
        a = self.net.labels(pores=[0, 1, 2], mode='difference')
        b = ['pore.back', 'pore.bottom', 'pore.right', 'pore.top']
        assert sorted(a) == sorted(b)

    def test_labels_pores_mode_none(self):
        a = self.net.labels(pores=[0, 1], mode='none')
        assert a[0] != a[1]

    def test_labels_pores_mode_foo(self):
        with pytest.raises(Exception):
            self.net.labels(pores=[0, 1], mode='foo')

    def test_labels_hidden_key(self):
        self.net['pore._foo'] = True
        assert 'pore._foo' not in self.net.__str__()
        assert 'pore._foo' in self.net.keys()

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

    def test_setitem_wrong_prefix(self):
        with pytest.raises(Exception):
            self.geo['pore2.test'] = 0

    def test_setitem_wrong_length(self):
        with pytest.raises(Exception):
            self.geo['pore.test'] = sp.ones((self.geo.Np+1))
        assert 'pore.test' not in self.geo.keys()

    def test_setitem_replace_all(self):
        array_len = sp.size(self.geo['pore.all'])
        self.geo['pore.all'] = sp.ones((self.geo.Np+1))
        assert sp.size(self.geo['pore.all']) == array_len

    def test_setitem_overwrite_into_all(self):
        pass
        # This test will fail as there is currently no way to prevent this
        # array_sum = sp.sum(self.geo['pore.all'])
        # self.geo['pore.all'][0] = False
        # assert sp.sum(self.geo['pore.all']) == array_sum

    def test_object_name_name_conflict(self):
        with pytest.raises(Exception):
            self.geo.name = self.net.name

    def test_object_name_array_conflict(self):
        with pytest.raises(Exception):
            self.geo.name = 'coords'
        Np = self.geo.Np
        Nt = self.geo.Nt
        assert self.geo.Np == Np
        assert self.geo.Nt == Nt

    def test_get_indices(self):
        temp = self.net.pop('pore.all')
        with pytest.raises(Exception):
            self.net._get_indices(element='pores', labels='blah')
        self.net.update({'pore.all': temp})

    def test_get_indices_wildcard(self):
        a = self.net._get_indices(element='pore', labels='ri*')
        assert sp.all(a == [6, 7, 8, 15, 16, 17, 24, 25, 26])
        b = self.net._get_indices(element='pore', labels='*ght')
        assert sp.all(a == b)

    def test_write_dict(self):
        self.net['pore.test_dict'] = {'test1': 1, 'test2': self.net.Ps}
        assert 'pore.test_dict.test1' in self.net.keys()
        assert self.net['pore.test_dict.test1'].shape == (self.net.Np, )
        assert 'pore.test_dict.test2' in self.net.keys()

    def test_pore_mapping(self):
        a = self.geo21['pore._id']
        b = self.geo22['pore._id']
        assert a.size == self.geo21.Np
        assert b.size == self.geo22.Np
        assert ~sp.any(sp.in1d(a, b))
        assert sp.all(self.net2.map_pores(a) == self.net2.pores(self.geo21.name))
        assert sp.all(self.net2.map_pores(b) == self.net2.pores(self.geo22.name))

    def test_throat_mapping(self):
        a = self.geo21['throat._id']
        assert a.size == self.geo21.Nt
        assert sp.all(self.net2.map_throats(a) == self.net2.throats(self.geo21.name))

    def test_map_pores_unfiltered(self):
        a = self.geo['pore._id']
        b = self.net.map_pores(a, filtered=False)
        assert sp.all(b.indices == self.net.pores(self.geo.name))
        assert b.mask.size == self.geo.Np

    def test_interleave_data_bool(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        Ps = net.pores('top')
        geom1 = op.geometry.GenericGeometry(network=net, pores=Ps)
        Ps = net.pores('bottom')
        geom2 = op.geometry.GenericGeometry(network=net, pores=Ps)
        # Ensure Falses return in missing places
        geom1['pore.blah'] = True
        assert sp.all(~geom2['pore.blah'])
        assert sp.sum(net['pore.blah']) == 4
        # Ensure all Trues returned now
        geom2['pore.blah'] = True
        assert sp.all(geom2['pore.blah'])
        assert sp.sum(net['pore.blah']) == 8

    def test_interleave_data_int(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        Ps = net.pores('top')
        geom1 = op.geometry.GenericGeometry(network=net, pores=Ps)
        Ps = net.pores('bottom')
        geom2 = op.geometry.GenericGeometry(network=net, pores=Ps)
        geom1['pore.blah'] = 1
        # Ensure ints are returned geom1
        assert 'int' in geom1['pore.blah'].dtype.name
        # Ensure nans are returned on geom2
        assert sp.all(sp.isnan(geom2['pore.blah']))
        # Ensure interleaved array is float with nans
        assert 'float' in net['pore.blah'].dtype.name
        # Ensure missing values are floats
        assert sp.sum(sp.isnan(net['pore.blah'])) == 4

    def test_interleave_data_float(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        Ps = net.pores('top')
        geom1 = op.geometry.GenericGeometry(network=net, pores=Ps)
        Ps = net.pores('bottom')
        geom2 = op.geometry.GenericGeometry(network=net, pores=Ps)
        geom1['pore.blah'] = 1.0
        # Ensure flaots are returned geom1
        assert 'float' in geom1['pore.blah'].dtype.name
        # Ensure nans are returned on geom2
        assert sp.all(sp.isnan(geom2['pore.blah']))
        # Ensure interleaved array is float with nans
        assert 'float' in net['pore.blah'].dtype.name
        # Ensure missing values are floats
        assert sp.sum(sp.isnan(net['pore.blah'])) == 4

    def test_interleave_data_object(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        Ps = net.pores('top')
        geom1 = op.geometry.GenericGeometry(network=net, pores=Ps)
        Ps = net.pores('bottom')
        geom2 = op.geometry.GenericGeometry(network=net, pores=Ps)
        geom1['pore.blah'] = [[1, 2], [1, 2, 3], [1, 2, 3, 4], [1]]
        assert 'object' in net['pore.blah'].dtype.name
        # Ensure missing elements are None
        assert sp.sum([item is None for item in net['pore.blah']]) == 4

    def test_interleave_data_key_error(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        Ps = net.pores('top')
        geom1 = op.geometry.GenericGeometry(network=net, pores=Ps)
        Ps = net.pores('bottom')
        geom2 = op.geometry.GenericGeometry(network=net, pores=Ps)
        with pytest.raises(KeyError):
            net['pore.blah']
        with pytest.raises(KeyError):
            geom1['pore.blah']
        with pytest.raises(KeyError):
            geom2['pore.blah']

    def test_interleave_data_float_missing_geometry(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        geom = op.geometry.GenericGeometry(network=net, pores=[0, 1, 2])
        geom['pore.blah'] = 1.0
        assert sp.any(sp.isnan(net['pore.blah']))

    def test_interleave_data_int_missing_geometry(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        geom = op.geometry.GenericGeometry(network=net, pores=[0, 1, 2])
        geom['pore.blah'] = 1
        assert sp.any(sp.isnan(net['pore.blah']))

    def test_interleave_data_bool_missing_geometry(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        geom = op.geometry.GenericGeometry(network=net, pores=[0, 1, 2])
        geom['pore.blah'] = True
        assert sp.sum(net['pore.blah']) == geom.Np

    def test_interpolate_data(self):
        a = self.geo.interpolate_data(propname='throat.diameter')
        assert a.size == self.geo.Np
        a = self.geo.interpolate_data(propname='pore.diameter')
        assert a.size == self.geo.Nt

#    def test_get_regenerate_on_demand(self):
#        self.geo.regenerate_models()
#        models = list(self.geo.models.keys())
#        assert len(set(self.geo.keys()).intersection(models)) == 2
#        for item in models:
#            del self.geo[item]
#        assert len(set(self.geo.keys()).intersection(models)) == 0
#        for item in models:
#            self.geo.models[item]['regen_mode'] = 'deferred'
#            arr = self.geo[item]
#            self.geo.models[item]['regen_mode'] = 'normal'
#        assert len(set(self.geo.keys()).intersection(models)) == 2
#
#    def test_get_no_matches(self):
#        self.geo.pop('pore.blah', None)
#        with pytest.raises(KeyError):
#            self.geo['pore.blah']


if __name__ == '__main__':

    t = BaseTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
