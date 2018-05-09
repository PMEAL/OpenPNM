import openpnm as op
import scipy as sp
import pytest


class TransientStokesFlowTest:

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
        mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
        self.phys.add_model(propname='throat.hydraulic_conductance',
                            model=mod,
                            viscosity='throat.viscosity',
                            regen_mode='normal')

    def test_transient_stokes_flow(self):
        alg = op.algorithms.TransientStokesFlow(network=self.net,
                                                phase=self.phase)
        alg.settings.update({'t_scheme': 'implicit', 't_step': 1,
                             't_output': 100, 't_tolerance': 1e-08})
        alg.set_IC(0)
        alg.set_value_BC(pores=self.net.pores('back'), values=1)
        alg.set_value_BC(pores=self.net.pores('front'), values=0)
        alg.run()
        x = [0., 0., 0.,
             0.22636, 0.37759, 0.43854,
             0.6402, 0.83987, 0.77547,
             1., 1., 1.]
        y = sp.around(alg[alg.settings['quantity']], decimals=5)
        assert sp.all(x == y)

    def teardown_class(self):
        ws = op.core.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = TransientStokesFlowTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
