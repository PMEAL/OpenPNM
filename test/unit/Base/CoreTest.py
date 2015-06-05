import OpenPNM
import scipy as sp

class CoreTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3,3,3], name='test_net')
        self.geo = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                    pores=self.net.Ps,
                                                    throats=self.net.Ts,
                                                    name='test_geo')
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
        assert sorted(a) == ['pore.diameter','pore.volume',
                             'throat.area','throat.diameter']

    def test_props_models(self):
        a = self.geo.props(mode='models')
        assert sorted(a) == ['pore.volume','throat.area']

    def test_props_constants(self):
        a = self.geo.props(mode='constants')
        assert sorted(a) == ['pore.diameter','throat.diameter']

    def test_props_pores_all(self):
        a = self.geo.props(element='pores')
        assert sorted(a) == ['pore.diameter','pore.volume']

    def test_props_pores_models(self):
        a = self.geo.props(element='pores', mode='models')
        assert sorted(a) == ['pore.volume']

    def test_props_pores_constants(self):
        a = self.geo.props(element='pores', mode='constants')
        assert sorted(a) == ['pore.diameter']

    def test_pores(self):
        a = self.net.pores()
        assert sp.all(a == sp.arange(0,self.net.Np))

    def test_pores_one_label(self):
        a = self.net.pores(labels='top')
        assert sp.all(a == [ 2,  5,  8, 11, 14, 17, 20, 23, 26])

    def test_pores_two_labels_union(self):
        a = self.net.pores(labels=['top','front'], mode='union')
        assert sp.all(a == [ 0,  1,  2,  3,  4,  5,  6,
                             7,  8, 11, 14, 17, 20, 23, 26])

    def test_pores_two_labels_intersection(self):
        a = self.net.pores(labels=['top','front'], mode='intersection')
        assert sp.all(a == [ 2, 5, 8])

    def test_pores_two_labels_not_intersection(self):
        a = self.net.pores(labels=['top','front'], mode='not_intersection')
        assert sp.all(a == [ 9, 10, 12, 13, 15, 16, 18, 19, 21, 22, 24, 25])

    def test_throats(self):
        a = self.net.throats()
        assert sp.all(a == sp.arange(0,self.net.Nt))

    def test_throats_one_label(self):
        a = self.net.throats(labels='label1')
        assert sp.all(a == [0, 1, 2, 3, 4, 5])

    def test_throats_two_labels_union(self):
        a = self.net.throats(labels=['label1','label2'], mode='union')
        assert sp.all(a == [0, 1, 2, 3, 4, 5, 6, 7, 8])

    def test_throats_two_labels_intersection(self):
        a = self.net.throats(labels=['label1','label2'], mode='intersection')
        assert sp.all(a == [ 3, 4, 5])

    def test_throats_two_labels_not_intersection(self):
        a = self.net.throats(labels=['label1','label2'], mode='not_intersection')
        assert sp.all(a == [0, 1, 2, 6, 7, 8])

    def test_tomask_pores(self):
        a = self.net.tomask(pores=self.net.pores('top'))
        assert sp.sum(a) == 9

    def test_tomask_throats(self):
        a = self.net.tomask(throats=self.net.throats('label1'))
        assert sp.sum(a) == 6

    def test_toindices_pores(self):
        mask = sp.zeros((self.net.Np,),dtype=bool)
        Ps = [0,3,6]
        mask[Ps] = True
        a = self.net.toindices(mask)
        assert sp.all(a == Ps)

    def test_toindices_throats(self):
        mask = sp.zeros((self.net.Nt,),dtype=bool)
        Ts = [0,3,6]
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
        a = self.net.num_pores(labels=['top','front'],mode='union')
        assert a == 15

    def test_num_pores_two_labels_intersection(self):
        a = self.net.num_pores(labels=['top','front'],mode='intersection')
        assert a == 3

    def test_num_pores_two_labels_notintersection(self):
        a = self.net.num_pores(labels=['top','front'],mode='not_intersection')
        assert a == 12

    def test_num_pores_two_labels_difference(self):
        a = self.net.num_pores(labels=['top','front'],mode='difference')
        assert a == 12

    def test_num_throats(self):
        a = self.net.num_throats()
        assert a == 54

    def test_num_throats_one_label(self):
        a = self.net.num_throats(labels='label1')
        assert a == 6

    def test_num_throats_two_labels_union(self):
        a = self.net.num_throats(labels=['label1','label2'],mode='union')
        assert a == 9

    def test_num_throats_two_labels_intersection(self):
        a = self.net.num_throats(labels=['label1','label2'],mode='intersection')
        assert a == 3

    def test_num_throats_two_labels_notintersection(self):
        a = self.net.num_throats(labels=['label1','label2'],mode='not_intersection')
        assert a == 6

    def test_num_throats_two_labels_difference(self):
        a = self.net.num_throats(labels=['label1','label2'],mode='difference')
        assert a == 45



