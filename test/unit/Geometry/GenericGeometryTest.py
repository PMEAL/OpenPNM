import OpenPNM
import scipy as sp
import pytest


class GenericGeometryTest:

    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        P = self.net.pores('top', mode='not')
        T = self.net.find_neighbor_throats(pores=P)
        self.geo = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                    pores=P,
                                                    throats=T)
        P = self.net.pores('top')
        T = self.net.find_neighbor_throats(pores=P, mode='intersection')
        self.geo2 = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                     pores=P,
                                                     throats=T)

    def test_init_without_network(self):
        standalone_geo = OpenPNM.Geometry.GenericGeometry()
        assert standalone_geo.Np == 0
        assert standalone_geo.Nt == 0

    def test_getitem_coords_from_net(self):
        a = self.geo['pore.coords']
        b = self.net['pore.coords'][self.net.pores(self.geo.name)]
        assert sp.all(a == b)

    def test_getitem_geo_conns_includes_nans(self):
        a = self.geo['throat.conns']
        assert sp.all(len(a) == self.geo.Nt)
        assert sp.amax(a) == (self.geo.Np - 1)
        b = [pore for pore in a[:, 1] if pore is not sp.nan]
        assert sp.size(b) > self.geo.Np - 25

    def test_getitem_geo2_conns_has_no_nans(self):
        # Because of the way geo2 was constructed it has not dangling throats
        # to another Geometry, so should have no nans in its conns
        a = self.geo2['throat.conns']
        assert sp.all(len(a) == self.geo2.Nt)
        assert sp.amax(a) == (self.geo2.Np - 1)
        b = sp.where(~(a[:, 1] < sp.inf))[0]
        assert sp.size(b) == 0

    def test_getitem_own_pores_by_name(self):
        a = self.geo.pores(self.geo.name)
        b = self.geo.Ps
        assert sp.all(a == b)

    def test_get_item_self_name(self):
        a = self.geo.get('pore.'+self.geo.name)
        assert a is None
        a = self.geo['pore.'+self.geo.name]
        assert sp.sum(a) == self.geo.Np

    def test_initialize_with_overlapping_locations(self):
        with pytest.raises(Exception):
            OpenPNM.Geometry.GenericGeometry(network=self.net, pores=[0])

    def test_add_and_remove_pores(self):
        Ps = self.geo.Pnet[:50]
        self.geo.set_locations(pores=Ps, mode='remove')
        assert self.geo.Np == 50
        self.geo.set_locations(pores=Ps, mode='add')
        assert self.geo.Np == 100

    def test_add_and_remove_pores_with_labels(self):
        self.geo['pore.label'] = False
        self.geo['pore.label'][40:80] = True
        assert self.geo.num_pores('label') == 40
        Ps = self.geo.Pnet[:50]
        self.geo.set_locations(pores=Ps, mode='remove')
        assert self.geo.num_pores('label') == 30
        self.geo.set_locations(pores=Ps, mode='add')
        assert self.geo.num_pores('label') == 30

    def test_plot_histogram(self):
        self.geo['pore.diameter'] = 1
        self.geo['throat.diameter'] = 1
        self.geo['throat.length'] = 1
        self.geo.plot_histograms()

    def test_clear(self):
        self.geo2.clear(mode='complete')
        assert len(self.geo2.props()) == 0
        assert len(self.geo2['pore.all']) == 0
        assert len(self.geo2['throat.all']) == 0
        assert self.net.num_pores(self.geo2.name) == 0
        assert self.net.num_throats(self.geo2.name) == 0
        # Repair geo2 for use in other tests?
        ctrl = self.net.controller
        ctrl.purge_object(self.geo2)
        P = self.net.pores('top')
        T = self.net.find_neighbor_throats(pores=P, mode='intersection')
        self.geo2 = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                     pores=P,
                                                     throats=T)
