import OpenPNM
import scipy as sp


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

    def test_props_all(self):
        a = self.geo.props()
        assert sorted(a) == ['pore.diameter', 'pore.volume',
                             'throat.area', 'throat.diameter']

    def test_props_models(self):
        a = self.geo.props(mode='models')
        assert sorted(a) == ['pore.volume', 'throat.area']

    def test_props_constants(self):
        a = self.geo.props(mode='constants')
        assert sorted(a) == ['pore.diameter', 'throat.diameter']

    def test_props_pores_all(self):
        a = self.geo.props(element='pores')
        assert sorted(a) == ['pore.diameter', 'pore.volume']

    def test_props_pores_models(self):
        a = self.geo.props(element='pores', mode='models')
        assert sorted(a) == ['pore.volume']

    def test_props_pores_constants(self):
        a = self.geo.props(element='pores', mode='constants')
        assert sorted(a) == ['pore.diameter']

    def test_labels(self):
        a = self.net.labels()
        assert 'pore.top' in a

    def test_labels_on_pores(self):
        a = self.net.labels(element='pores')
        assert sorted(a) == ['pore.all', 'pore.back', 'pore.bottom',
                             'pore.front', 'pore.internal', 'pore.left',
                             'pore.right', 'pore.'+self.geo.name, 'pore.top']

    def test_labels_on_throats(self):
        a = self.net.labels(element='throats')
        assert sorted(a) == ['throat.all', 'throat.'+self.geo.name]

    def test_labels_on_foo(self):
        a = self.net.labels(element='foo')
        assert a is None

    def test_labels_on_all_pores(self):
        a = self.net.labels(pores='all')
        assert sorted(a) == ['pore.all', 'pore.back', 'pore.bottom',
                             'pore.front', 'pore.internal', 'pore.left',
                             'pore.right', 'pore.'+self.geo.name, 'pore.top']

    def test_labels_on_all_throats(self):
        a = self.net.labels(throats='all')
        assert sorted(a) == ['throat.all', 'throat.'+self.geo.name]

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

    def test_clear(self):
        Np = self.geo.Np
        Nt = self.geo.Nt
        self.geo.clear()
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
        assert a == b
