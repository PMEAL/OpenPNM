import openpnm as op
import numpy as np


class BoundaryTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        Ps_int = self.net.pores(labels=['top', 'bottom'], mode='not')
        Ps_boun = self.net.pores(labels=['top', 'bottom'], mode='union')
        Pb_mask = np.random.random(len(Ps_boun)) < 0.5
        Ts_int = self.net.throats(labels=['internal'])
        Ts_boun = self.net.throats(labels=['internal'], mode='not')
        Tb_mask = np.random.random(len(Ts_boun)) < 0.5
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=Ps_int,
                                            throats=Ts_int)
        self.boun1 = op.geometry.Boundary(network=self.net,
                                          shape='spheres',
                                          pores=Ps_boun[Pb_mask],
                                          throats=Ts_boun[Tb_mask])
        self.boun2 = op.geometry.Boundary(network=self.net,
                                          shape='cubes',
                                          pores=Ps_boun[~Pb_mask],
                                          throats=Ts_boun[~Tb_mask])
        self.geo.regenerate_models()

    def teardown_class(self):
        mgr = op.Workspace()
        mgr.clear()

    def test_plot_histogram(self):
        for obj in [self.boun1, self.boun2]:
            obj.show_hist()
            obj.show_hist(props=['pore.diameter', 'pore.volume',
                                 'throat.length'])
            obj.show_hist(props=['pore.diameter', 'pore.volume',
                                 'throat.length', 'throat.diameter',
                                 'pore.seed'])


if __name__ == '__main__':

    t = BoundaryTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
