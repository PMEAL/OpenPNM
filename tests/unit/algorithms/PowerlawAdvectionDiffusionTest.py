import openpnm as op
import scipy as sp
import pytest


class PowerlawAdvectionDiffusionTest:

    def setup_class(self):
        sp.random.seed(0)
        self.net = op.network.Cubic(shape=[4, 3, 1], spacing=1e-4)
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)

        self.phase = op.phases.Water(network=self.net)
        self.phase['throat.viscosity'] = self.phase['pore.viscosity'][0]

        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        mod1 = op.models.physics.hydraulic_conductance.hagen_poiseuille
        self.phys.add_model(propname='throat.hydraulic_conductance',
                            model=mod1, viscosity='throat.viscosity',
                            regen_mode='normal')
        mod2 = op.models.physics.diffusive_conductance.ordinary_diffusion
        self.phys.add_model(propname='throat.diffusive_conductance',
                            model=mod2, regen_mode='normal')

    def test_powerlaw_advection_diffusion_diffusion(self):
        alg1 = op.algorithms.StokesFlow(network=self.net, phase=self.phase)
        alg1.set_value_BC(pores=self.net.pores('back'), values=10)
        alg1.set_value_BC(pores=self.net.pores('front'), values=0)
        alg1.run()
        self.phase[alg1.settings['quantity']] = alg1[alg1.settings['quantity']]

        alg2 = op.algorithms.AdvectionDiffusion(network=self.net,
                                                phase=self.phase)
        alg2.settings.update({'s_scheme': 'powerlaw'})
        alg2.set_value_BC(pores=self.net.pores('back'), values=2)
        alg2.set_value_BC(pores=self.net.pores('front'), values=0)
        alg2.run()
        x = [0., 0., 0.,
             0.78223, 0.97971, 1.06055,
             1.50462, 1.73478, 1.60123,
             2., 2., 2.]
        y = sp.around(alg2[alg2.settings['quantity']], decimals=5)
        assert sp.all(x == y)

    def teardown_class(self):
        ws = op.core.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = PowerlawAdvectionDiffusionTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
