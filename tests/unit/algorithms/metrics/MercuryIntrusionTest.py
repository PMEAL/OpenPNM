import numpy as np
import openpnm as op
mgr = op.Workspace()


class MercuryIntrusionTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[15, 15, 15], spacing=0.0005)
        self.geo = op.geometry._StickAndBall(network=self.net,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)

    def test_run(self):
        mip = op.algorithms.metrics.MercuryIntrusion(network=self.net)
        mip.run()
        assert mip['pore.invasion_sequence'].max() > 2
        assert mip['throat.invasion_sequence'].max() > 2

    def test_pcsnwp_data(self):
        mip = op.algorithms.metrics.MercuryIntrusion(network=self.net)
        mip.run()
        mip.pc_data = [1e1, 1e3, 1e5, 1e6]
        mip.snwp_data = [0.0, 0.1, 0.5, 0.9]
        assert np.all(mip.pc_data == [1e1, 1e3, 1e5, 1e6])
        assert np.all(mip.snwp_data == [0.0, 0.1, 0.5, 0.9])

    def test_plot_data(self):
        mip = op.algorithms.metrics.MercuryIntrusion(network=self.net)
        mip.run()
        mip.pc_data = [1e1, 1e3, 1e5, 1e6]
        mip.snwp_data = [0.0, 0.1, 0.5, 0.9]
        mip.plot_intrusion_curve()


if __name__ == '__main__':

    t = MercuryIntrusionTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
