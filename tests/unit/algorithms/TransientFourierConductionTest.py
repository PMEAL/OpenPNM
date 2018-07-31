import openpnm as op
import scipy as sp
import pytest


class TransientFourierConductionTest:

    def setup_class(self):
        sp.random.seed(0)
        self.net = op.network.Cubic(shape=[4, 3, 1], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)

        self.phase = op.phases.GenericPhase(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['throat.thermal_conductance'] = 1e-15
        self.geo['pore.volume'] = 1e-27

    def test_transient_fickian_diffusion(self):
        alg = op.algorithms.TransientFourierConduction(network=self.net,
                                                       phase=self.phase)
        alg.setup(quantity='pore.temperature',
                  conductance='throat.thermal_conductance',
                  t_initial=0, t_final=1000, t_step=1, t_output=100,
                  t_tolerance=1e-08, t_scheme='implicit')
        alg.set_IC(0)
        alg.set_value_BC(pores=self.net.pores('back'), values=1)
        alg.set_value_BC(pores=self.net.pores('front'), values=0)
        alg.run()
        x = [0., 0., 0.,
             0.33333, 0.33333, 0.33333,
             0.66667, 0.66667, 0.66667,
             1., 1., 1.]
        y = sp.around(alg[alg.settings['quantity']], decimals=5)
        assert sp.all(x == y)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = TransientFourierConductionTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
