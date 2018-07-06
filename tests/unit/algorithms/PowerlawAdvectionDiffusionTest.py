import openpnm as op
import scipy as sp
import pytest


class PowerlawAdvectionDiffusionTest:

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
        self.phys['throat.diffusive_conductance'] = 1e-15
        self.phys['throat.hydraulic_conductance'] = 1e-15

    def test_powerlaw_advection_diffusion_diffusion(self):
        sf = op.algorithms.StokesFlow(network=self.net, phase=self.phase)
        sf.set_value_BC(pores=self.net.pores('back'), values=1)
        sf.set_value_BC(pores=self.net.pores('front'), values=0)
        sf.run()
        self.phase[sf.settings['quantity']] = sf[sf.settings['quantity']]

        ad = op.algorithms.AdvectionDiffusion(network=self.net, phase=self.phase)
        ad.settings.update({'s_scheme': 'powerlaw'})
        ad.set_value_BC(pores=self.net.pores('back'), values=2)
        ad.set_value_BC(pores=self.net.pores('front'), values=0)
        ad.run()

        x = [0., 0., 0.,
             0.89653, 0.89653, 0.89653,
             1.53924, 1.53924, 1.53924,
             2., 2., 2.]
        y = sp.around(ad[ad.settings['quantity']], decimals=5)
        assert sp.all(x == y)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = PowerlawAdvectionDiffusionTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
