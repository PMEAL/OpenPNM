import OpenPNM


class Toray090Test:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.geo = OpenPNM.Geometry.Toray090(network=self.net,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)

    def test_props_exist(self):
        props = ['pore.seed',
                 'pore.diameter',
                 'pore.volume',
                 'pore.area',
                 'throat.length',
                 'throat.seed',
                 'throat.diameter',
                 'throat.volume',
                 'throat.area']
        assert set(props).issubset(set(self.geo.keys()))
