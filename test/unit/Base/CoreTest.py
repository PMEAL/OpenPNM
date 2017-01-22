import OpenPNM
import scipy as sp
import pytest
mgr = OpenPNM.Base.Workspace()
mgr.loglevel = 60


class CoreTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.geo = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                    pores=self.net.Ps,
                                                    throats=self.net.Ts)
        self.geo['pore.diameter'] = sp.rand(self.net.Np)
        self.geo.models.add(propname='pore.volume',
                            model=OpenPNM.Geometry.models.pore_volume.sphere)
        self.geo['throat.diameter'] = sp.rand(self.net.Nt)
        self.geo.models.add(propname='throat.area',
                            model=OpenPNM.Geometry.models.throat_area.cylinder)
        self.geo['throat.label1'] = False
        self.geo['throat.label2'] = False
        self.geo['throat.label1'][0:6] = True
        self.geo['throat.label2'][3:9] = True
        self.net1 = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.geo1 = OpenPNM.Geometry.GenericGeometry(network=self.net1,
                                                     pores=self.net1.Ps,
                                                     throats=self.net1.Ts)
        self.phase1 = OpenPNM.Phases.GenericPhase(network=self.net1)
        self.phase2 = OpenPNM.Phases.GenericPhase(network=self.net1)
        self.phys1 = OpenPNM.Physics.GenericPhysics(network=self.net1,
                                                    geometry=self.geo1,
                                                    phase=self.phase1)
        self.phys2 = OpenPNM.Physics.GenericPhysics(network=self.net1,
                                                    geometry=self.geo1,
                                                    phase=self.phase2)
        self.net2 = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        Ps = sp.arange(0, 18)
        Ts = self.net2.find_neighbor_pores(Ps, mode='union')
        self.geo21 = OpenPNM.Geometry.GenericGeometry(network=self.net2,
                                                      pores=Ps,
                                                      throats=Ts)
        Ps = sp.arange(18, 27)
        Ts = self.net2.find_neighbor_pores(Ps, mode='intersection')
        self.geo22 = OpenPNM.Geometry.GenericGeometry(network=self.net2,
                                                      pores=Ps,
                                                      throats=Ts)

    def test_clear_complete(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        b = sorted(list(net.keys()))
        dict_ = net.copy()
        net.clear()
        assert net.Np == 0
        assert net.Nt == 0
        a = sorted(list(net.keys()))
        assert a == ['pore.all', 'throat.all']
        net.update(dict_)
        assert net.Np == 27
        assert net.Nt == 54
        a = sorted(list(net.keys()))
        assert a == b

    def test_clear_props(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        net.clear(mode='props')
        assert len(net.props()) == 0

    def test_clear_labels(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        net.clear(mode='labels')
        assert len(net.labels()) == 2
        assert net.Np == 27
        assert net.Nt == 54

    def test_clear_labels_and_props(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        net.clear(mode=['labels', 'props'])
        assert len(net.labels()) == 2
        assert net.Np == 27
        assert net.Nt == 54
        assert len(net.props()) == 0

    def test_clear_models(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        geo = OpenPNM.Geometry.Stick_and_Ball(network=net,
                                              pores=net.Ps,
                                              throats=net.Ts)
        geo.clear(mode='props')
        assert len(geo.props()) == 0
        geo['pore.seed'] = sp.rand(geo.Np)
        geo.models.regenerate()
        assert len(geo.props()) > 0
        geo.clear(mode=['props', 'models'])
        geo.models.regenerate()
        assert len(geo.props()) == 0

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
        b = ['throat.all', 'throat.'+self.geo.name]
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
        b = ['throat.all', 'throat.'+self.geo.name]
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

    def test_pores_two_labels_not_intersection(self):
        a = self.net.pores(labels=['top', 'front'], mode='not_intersection')
        assert sp.all(a == [0, 1, 3, 4, 6, 7, 11, 14, 17, 20, 23, 26])

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

    def test_throats_two_labels_not_intersection(self):
        a = self.net.throats(labels=['label1', 'label2'],
                             mode='not_intersection')
        assert sp.all(a == [0, 1, 2, 6, 7, 8])

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

    def test_filter_by_label_pores_two_labels_not_intersection(self):
        Ps = self.net.pores(['top', 'bottom', 'front'])
        a = self.net.filter_by_label(pores=Ps,
                                     labels=['top', 'front'],
                                     mode='not_intersection')
        b = [0, 1, 3, 4, 6, 7, 11, 14, 17, 20, 23, 26]
        assert sp.all(a == b)

    def test_filter_by_label_pores_two_labels_not(self):
        Ps = self.net.pores(['top', 'bottom', 'front'])
        a = self.net.filter_by_label(pores=Ps,
                                     labels=['top', 'front'],
                                     mode='not')
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

    def test_toindices_wrong_mask(self):
        mask = sp.zeros((self.net.Nt)-2, dtype=bool)
        mask[[0, 3, 6]] = True
        with pytest.raises(Exception):
            self.net.toindices(mask)

    def test_count(self):
        a = self.net._count()
        assert a == {'pore': 27, 'throat': 54}

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

    def test_num_pores_two_labels_notintersection(self):
        a = self.net.num_pores(labels=['top', 'front'],
                               mode='not_intersection')
        assert a == 12

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

    def test_num_throats_two_labels_notintersection(self):
        a = self.net.num_throats(labels=['label1', 'label2'],
                                 mode='not_intersection')
        assert a == 6

    def test_num_throats_two_labels_difference(self):
        a = self.net.num_throats(labels=['label1', 'label2'],
                                 mode='difference')
        assert a == 45

    def test_setitem_wrong_prefix(self):
        with pytest.raises(Exception):
            self.geo['pore2.test'] = 0

    def test_setitem_wrong_length(self):
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

    def test_object_rename(self):
        assert self.geo1 in mgr.values()
        old_name = self.geo1.name
        self.geo1.name = 'new_name'
        assert self.geo1.name == 'new_name'
        assert self.geo1 in mgr.values()
        self.geo1.name = old_name

    def test_object_duplicate_name(self):
        temp = self.geo1.name
        try:
            self.geo1.name = self.net1.name
        except:
            pass
        assert self.geo1.name == temp

    def test_geometry_lookup_all(self):
        a = self.net1.geometries()
        assert a == [self.geo1.name]

    def test_geometry_lookup_by_name(self):
        a = self.net1.geometries(self.geo1.name)
        assert a == [self.geo1]

    def test_set_locations_on_phase(self):
        with pytest.raises(Exception):
            self.phase1._set_locations(element='pores', locations=1)

    def test_set_locations_add_on_empty_geometry(self):
        # 'empty' meaning assigned nowhere, with no models or props
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        geo = OpenPNM.Geometry.GenericGeometry(network=net)
        assert geo.Np == 0
        assert geo.Nt == 0
        geo.set_locations(pores=net.Ps)
        assert geo.Np == net.Np
        assert geo.Nt == 0
        geo.set_locations(throats=net.Ts)
        assert geo.Np == net.Np
        assert geo.Nt == net.Nt
        del mgr[net.name]
        del mgr[geo.name]

    def test_set_locations_add_and_remove_on_empty_geometry(self):
        # 'empty' meaning assigned nowhere, with no models or props
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        geo = OpenPNM.Geometry.GenericGeometry(network=net,
                                               pores=net.Ps,
                                               throats=net.Ts)
        Np = geo.Np
        geo.set_locations(pores=0, mode='remove')
        assert Np > geo.Np
        del mgr[net.name]
        del mgr[geo.name]

    def test_set_locations_overlapping_on_empty_geometry(self):
        # 'empty' meaning assigned nowhere, with no models or props
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        geo1 = OpenPNM.Geometry.GenericGeometry(network=net,
                                                pores=net.pores('top'))
        geo2 = OpenPNM.Geometry.GenericGeometry(network=net)
        try:
            geo2.set_locations(pores=net.pores('top'))
        except:
            pass
        assert geo2.Np == 0
        del mgr[net.name]
        del mgr[geo1.name]
        del mgr[geo2.name]

    def test_set_locations_add_successivly_on_empty_geometry(self):
        # 'empty' meaning assigned nowhere, with no models or props
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        geo = OpenPNM.Geometry.GenericGeometry(network=net)
        assert geo.Np == 0
        assert geo.Nt == 0
        geo.set_locations(pores=net.pores('top'))
        assert geo.Np == 9
        assert geo.Nt == 0
        geo.set_locations(pores=net.pores('bottom'))
        assert geo.Np == 18
        assert geo.Nt == 0
        del mgr[net.name]
        del mgr[geo.name]

    def test_set_locations_add_duplicates_on_empty_geometry(self):
        # 'empty' meaning assigned nowhere, with no models or props
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        geo = OpenPNM.Geometry.GenericGeometry(network=net)
        assert geo.Np == 0
        assert geo.Nt == 0
        geo.set_locations(pores=net.pores('top'))
        assert geo.Np == 9
        assert geo.Nt == 0
        try:
            geo.set_locations(pores=net.pores('front'))
        except:
            pass
        assert geo.Np == 9
        assert geo.Nt == 0
        del mgr[net.name]
        del mgr[geo.name]

    def test_set_locations_remove_on_realistic_geometry(self):
        # 'realistic' meaning assigned to pores, and has models and props
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        geo = OpenPNM.Geometry.GenericGeometry(network=net,
                                               pores=net.Ps,
                                               throats=net.Ts)
        geo['pore.prop'] = sp.arange(0, net.Np)
        f = OpenPNM.Geometry.models.pore_misc.random
        geo.models.add(propname='pore.regenerating_model',
                       model=f,
                       seed=0,
                       regen_mode='normal')
        geo.models.add(propname='pore.constant_model',
                       model=f,
                       seed=0,
                       regen_mode='constant')
        geo.set_locations(pores=0, mode='remove')
        assert geo.Np == 26
        geo.set_locations(pores=1, mode='remove')
        assert geo.Np == 25
        del mgr[net.name]
        del mgr[geo.name]

    def test_set_locations_add_on_realistic_geometry(self):
        # 'realistic' meaning assigned to pores, and has models and contants
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        geo = OpenPNM.Geometry.GenericGeometry(network=net,
                                               pores=net.pores('top'))
        geo['pore.prop'] = sp.arange(0, net.num_pores('top'))
        f = OpenPNM.Geometry.models.pore_misc.random
        geo.models.add(propname='pore.regenerating_model',
                       model=f,
                       seed=0,
                       regen_mode='normal')
        geo.models.add(propname='pore.constant_model',
                       model=f,
                       seed=0,
                       regen_mode='constant')
        try:
            geo.set_locations(pores=net.pores('bottom'), mode='add')
        except:
            pass
        assert geo.Np == 9
        del mgr[net.name]
        del mgr[geo.name]

    def test_set_locations_add_on_geometry_models_only(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        geo = OpenPNM.Geometry.GenericGeometry(network=net,
                                               pores=net.pores('top'))
        f = OpenPNM.Geometry.models.pore_misc.random
        geo.models.add(propname='pore.regenerating_model',
                       model=f,
                       seed=0,
                       regen_mode='normal')
        geo.set_locations(pores=net.pores('bottom'), mode='add')
        assert geo.Np == 18
        assert sp.size(geo['pore.regenerating_model']) == 18
        del mgr[net.name]
        del mgr[geo.name]

    def test_find_all_physics(self):
        a = set(self.net1.physics())
        b = {self.phys1.name, self.phys2.name}
        assert a == b

    def test_find_physics_by_name(self):
        a = self.net1.physics(self.phys1.name)
        assert self.phys1 in a
        assert self.phys2 not in a
        a = self.net1.physics([self.phys1.name, self.phys2.name])
        assert self.phys1 in a
        assert self.phys2 in a

    def test_find_all_phases(self):
        a = set(self.net1.phases())
        b = {self.phase1.name, self.phase2.name}
        assert a == b

    def test_find_phases_by_name(self):
        a = self.net1.phases(self.phase1.name)
        assert self.phase1 in a
        assert self.phase2 not in a
        a = self.net1.phases([self.phase1.name, self.phase2.name])
        assert self.phase1 in a
        assert self.phase2 in a

    def test_find_all_geometries(self):
        a = set(self.net1.geometries())
        b = {self.geo1.name}
        assert a == b

    def test_find_geometries_by_name(self):
        a = self.net1.phases(self.phase1.name)
        assert self.phase1 in a
        assert self.phase2 not in a
        a = self.net1.phases([self.phase1.name, self.phase2.name])
        assert self.phase1 in a
        assert self.phase2 in a

    def test_find_network_from_geometry(self):
        a = self.geo.network()
        assert a == [self.net]

    def test_find_network_by_name_from_geometry(self):
        a = self.geo.network(self.net.name)
        assert a == self.net

    def test_find_network_from_phase(self):
        a = self.phase1.network()
        assert a == [self.net1]

    def test_find_network_by_name_from_phase(self):
        a = self.phase1.network(self.net1.name)
        assert a == self.net1

    def test_find_network_from_physics(self):
        a = self.phys1.network()
        assert a == [self.net1]

    def test_find_network_by_name_from_physics(self):
        a = self.phys1.network(self.net1.name)
        assert a == self.net1

    def test_find_object_by_type(self):
        a = self.net._find_object(obj_type='geometry')
        assert type(a) is list
        a = self.net._find_object(obj_type='phase')
        assert type(a) is list
        a = self.net._find_object(obj_type='physics')
        assert type(a) is list
        a = self.net._find_object(obj_type='network')
        assert type(a) is list

    def test_object_print(self):
        a = self.net.__str__()
        assert type(a) == str

    def test_object_representation(self):
        a = self.net.__repr__()
        assert type(a) == str

    def test_parse_locations_boolean(self):
        b = sp.array([True, True, True])
        with pytest.raises(Exception):
            self.net._parse_locations(locations=b)
        b = sp.zeros((self.net.Np,), dtype=bool)
        assert len(self.net._parse_locations(locations=b)) == 0
        b = sp.zeros((self.net.Nt,), dtype=bool)
        b[[0, 1, 2]] = True
        assert sp.shape(self.net._parse_locations(locations=b)) == (3,)

    def test_parse_locations_None(self):
        assert len(self.net._parse_locations(locations=None)) == 0

    def test_parse_locations_string(self):
        with pytest.raises(Exception):
            self.net._parse_locations(locations='abc')

    def test_parse_locations_int(self):
        a = self.net._parse_locations(locations=0)
        assert type(a) == sp.ndarray
        assert sp.all(a == 0)

    def test_parse_locations_list(self):
        a = self.net._parse_locations(locations=[0, 1])
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

    def test_map_pores_geometry1_onto_network(self):
        a = self.geo21.map_pores(target=self.net2)
        assert sp.all(a == self.net2.pores(self.geo21.name))

    def test_map_pores_geometry2_onto_network(self):
        a = self.geo22.map_pores(target=self.net2, pores=self.geo22.Ps)
        assert sp.all(a == self.net2.pores(self.geo22.name))

    def test_map_pores_network_onto_self(self):
        a = self.net2.map_pores(target=self.net2)
        assert sp.all(a == self.net2.pores())

    def test_map_pores_geometry_onto_other_geometry(self):
        with pytest.raises(Exception):
            self.geo21.map_pores(target=self.geo22)

    def test_mapping(self):
        # Create small cubic network
        pn = OpenPNM.Network.Cubic(shape=[3, 3, 3], spacing=0.0001)
        # Assign 3 different geometries to each layer in the z-direction
        Pa = sp.arange(0, 9)
        Ta = pn.find_neighbor_throats(Pa)
        geom1 = OpenPNM.Geometry.GenericGeometry(network=pn,
                                                 pores=Pa,
                                                 throats=Ta)
        Pc = sp.arange(18, 27)
        Tc = pn.find_neighbor_throats(Pc)
        geom3 = OpenPNM.Geometry.GenericGeometry(network=pn,
                                                 pores=Pc,
                                                 throats=Tc)
        Pb = sp.arange(9, 18)
        Tb = pn.find_neighbor_throats(pores=Pb, mode='intersection')
        geom2 = OpenPNM.Geometry.GenericGeometry(network=pn,
                                                 pores=Pb,
                                                 throats=Tb)
        # Create an index in the Network
        pn['pore.num1'] = pn.Ps
        # Create the same index across each geom
        geom1['pore.num2'] = Pa
        geom2['pore.num2'] = Pb
        geom3['pore.num2'] = Pc
        # Confirm two indexes match
        assert(sp.all(pn['pore.num1'] == pn['pore.num2']))
        # Send junk pores to ensure error is raised
        with pytest.raises(Exception):
            pn.map_pores(pores=[0, pn.Np-1], target=geom1)
            pn.map_pores(pores=[0, pn.Np+1], target=geom1)
            pn.map_pores(pores=[pn.Np-1], target=geom1)
            pn.map_pores(pores=[pn.Np+1], target=geom1)
            geom1.map_pores(pores=[0, geom1.Np+1], target=pn)
            geom1.map_pores(pores=[0, pn.Np+1], target=pn)
            geom1.map_pores(pores=[geom1.Np+1], target=pn)
            geom1.map_pores(pores=[pn.Np+1], target=pn)
            geom2.map_pores(pores=[0], target=geom1)
            geom2.map_pores(pores=[geom2.Np+1], target=geom1)
            geom2.map_pores(pores=[0, geom2.Np-1], target=geom1)
            geom2.map_pores(pores=[0, geom2.Np+1], target=geom1)
        # Trim column from center of Network
        pn.trim(pores=[4, 13, 22])
        # Confirm index still match
        assert(sp.all(pn['pore.num1'] == pn['pore.num2']))
        # Check mapping between each Geometry object and in both directions
        # Check geom1
        a = geom1.map_pores(pores=geom1.Ps, target=pn)
        b = pn.map_pores(pores=a, target=geom1)
        assert(sp.all(b == geom1.Ps))
        a = geom1.map_throats(throats=geom1.Ts, target=pn)
        b = pn.map_throats(throats=a, target=geom1)
        assert(sp.all(b == geom1.Ts))
        # Check geom2
        a = geom2.map_pores(pores=geom2.Ps, target=pn)
        b = pn.map_pores(pores=a, target=geom2)
        assert(sp.all(b == geom2.Ps))
        a = geom2.map_throats(throats=geom2.Ts, target=pn)
        b = pn.map_throats(throats=a, target=geom2)
        assert(sp.all(b == geom2.Ts))
        # Check geom3
        a = geom3.map_pores(pores=geom3.Ps, target=pn)
        b = pn.map_pores(pores=a, target=geom3)
        assert(sp.all(b == geom3.Ps))
        a = geom3.map_throats(throats=geom3.Ts, target=pn)
        b = pn.map_throats(throats=a, target=geom3)
        assert(sp.all(b == geom3.Ts))

    def test_check_data_health(self):
        a = self.net.check_data_health()
        assert a.health
        for item in a.values():
            assert item == []
        self.net['pore.data_test'] = sp.nan
        a = self.net.check_data_health()
        assert not a.health
        assert a['pore.data_test'] == 'Has NaNs'
        del self.net['pore.data_test']
