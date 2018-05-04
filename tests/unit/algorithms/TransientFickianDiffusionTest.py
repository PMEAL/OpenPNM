import openpnm as op
import scipy as sp
import pytest


class TransientFickianDiffusionTest:

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
        mod = op.models.physics.diffusive_conductance.bulk_diffusion
        self.phys.add_model(propname='throat.diffusive_conductance',
                            model=mod, regen_mode='normal')

    def test_transient_fickian_diffusion(self):
        alg = op.algorithms.TransientFickianDiffusion(network=self.net,
                                                      phase=self.phase)
        alg.settings.update({'t_scheme': 'implicit', 't_step': 1,
                             't_output': 100, 't_tolerance': 1e-08})
        alg.set_IC(0)
        alg.set_dirichlet_BC(pores=self.net.pores('back'), values=1)
        alg.set_dirichlet_BC(pores=self.net.pores('front'), values=0)
        alg.run()
        x = [0., 0., 0.,
             0.28421, 0.35742, 0.3742,
             0.65712, 0.76651, 0.70061,
             1., 1., 1.]
        y = sp.around(alg[alg.settings['quantity']], decimals=5)
        assert sp.all(x == y)

    def teardown_class(self):
        ws = op.core.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = TransientFickianDiffusionTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
