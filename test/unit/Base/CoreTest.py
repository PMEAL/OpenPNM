import OpenPNM
import scipy as sp
ctrl = OpenPNM.Base.Controller()
ctrl.loglevel = 60


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
        a = self.net.labels(element='foo')
        assert a is None

    def test_labels_on_all_pores(self):
        a = self.net.labels(pores='all')
        b = ['pore.all', 'pore.back', 'pore.bottom', 'pore.front',
             'pore.internal', 'pore.left', 'pore.right', 'pore.'+self.geo.name,
             'pore.top']
        assert sorted(a) == sorted(b)

    def test_labels_on_all_throats(self):
        a = self.net.labels(throats='all')
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
        a = self.net.labels(pores=[0, 1], mode='foo')
        assert a is None

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
        a = self.net.filter_by_label(pores=Ps)
        assert sp.all(a == Ps.tolist())

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
        a = None
        try:
            a = self.net.toindices(mask)
        except:
            pass
        assert a is None

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
        self.geo['pore2.test'] = 0
        assert 'pore2.test' not in self.geo.keys()

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
        flag = False
        try:
            self.geo.name = self.net.name
        except:
            flag = True
        assert flag

    def test_object_name_array_conflict(self):
        flag = False
        try:
            self.geo.name = 'coords'
        except:
            flag = True
        assert flag

        Np = self.geo.Np
        Nt = self.geo.Nt
        assert self.geo.Np == Np
        assert self.geo.Nt == Nt

    def test_get_indices(self):
        flag = False
        temp = self.net.pop('pore.all')
        try:
            self.net._get_indices(element='pores', labels='blah')
        except:
            flag = True
        assert flag
        self.net.update({'pore.all': temp})

    def test_get_indices_wildcard(self):
        a = self.net._get_indices(element='pore', labels='ri*')
        assert sp.all(a == [6, 7, 8, 15, 16, 17, 24, 25, 26])
        b = self.net._get_indices(element='pore', labels='*ght')
        assert sp.all(a == b)

    def test_object_rename(self):
        assert self.geo1 in ctrl.values()
        self.geo1.name = 'new_name'
        assert self.geo1.name == 'new_name'
        assert self.geo1 in ctrl.values()

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
        del ctrl[net.name]
        del ctrl[geo.name]

    def test_set_locations_add_and_remove_on_empty_geometry(self):
        # 'empty' meaning assigned nowhere, with no models or props
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        geo = OpenPNM.Geometry.GenericGeometry(network=net,
                                               pores=net.Ps,
                                               throats=net.Ts)
        Np = geo.Np
        geo.set_locations(pores=0, mode='remove')
        assert Np > geo.Np
        del ctrl[net.name]
        del ctrl[geo.name]

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
        del ctrl[net.name]
        del ctrl[geo1.name]
        del ctrl[geo2.name]

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
        del ctrl[net.name]
        del ctrl[geo.name]

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
        del ctrl[net.name]
        del ctrl[geo.name]

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
        del ctrl[net.name]
        del ctrl[geo.name]

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
        del ctrl[net.name]
        del ctrl[geo.name]

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
        del ctrl[net.name]
        del ctrl[geo.name]

if __name__ == '__main__':
    a = CoreTest()
    a.setup_class()
    b = a.__class__.__dict__
    for item in b:
        if item.split('_')[0] == 'test':
            print('-'*79, '\n', item)
            b[item](self=a)
