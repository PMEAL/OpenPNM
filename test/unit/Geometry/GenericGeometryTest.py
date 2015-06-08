import OpenPNM
import scipy as sp


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
        b = sp.where(~(a[:, 1] < sp.inf))[0]
        assert sp.size(b) > 0

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
