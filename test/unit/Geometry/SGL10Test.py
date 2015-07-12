import OpenPNM


class SGL10Test:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.geo = OpenPNM.Geometry.SGL10(network=self.net,
                                          pores=self.net.Ps,
                                          throats=self.net.Ts)
