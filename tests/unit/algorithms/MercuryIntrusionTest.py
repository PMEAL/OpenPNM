import openpnm as op
import scipy as sp
mgr = op.Workspace()


class MercuryIntrusionTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[15, 15, 15], spacing=0.0005)
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)

    def test_run(self):
        mip = op.algorithms.MercuryIntrusion(network=self.net)
        assert mip['pore.invasion_sequence'].max() > 20
        assert mip['throat.invasion_sequence'].max() > 20

    def test_plot_pcsnwp_data(self):
        mip = op.algorithms.MercuryIntrusion(network=self.net)
        mip.pc_data = [1e1, 1e3, 1e5, 1e6]
        mip.snwp_data = [0.0, 0.1, 0.5, 0.9]
        assert sp.all(mip.pc_data == [1e1, 1e3, 1e5, 1e6])
        assert sp.all(mip.snwp_data == [0.0, 0.1, 0.5, 0.9])


if __name__ == '__main__':

    t = MercuryIntrusionTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
