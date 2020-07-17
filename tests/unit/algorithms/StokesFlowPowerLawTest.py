import openpnm as op
import numpy as np
from numpy.testing import assert_allclose


class StokesFlowPowerLawTest:

    def setup_class(self):
        np.random.seed(7)
        self.net = op.network.Cubic(shape=[5, 5, 1])

        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)

        self.phase = op.phases.Water(network=self.net)
        self.phase['pore.consistency'] = 4.2e-2  # Pa.s^n
        self.phase['pore.flow_index'] = 0.52
        self.phase['pore.viscosity_min'] = 0.001
        self.phase['pore.viscosity_max'] = 100
        self.phase['pore.viscosity'] = 0.0011

        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.flow = op.algorithms.NonNewtonianStokesFlow(network=self.net,
                                                         phase=self.phase)
        settings = {'cache_A': False, 'cache_b': False,
                    'solver_type': 'spsolve', 'relaxation_quantity': 0.7,
                    'max_iter': 200, 'solver_tol': 1e-6,
                    'iterative_props': (
                        'throat.nonNewtonian_hydraulic_conductance')}
        self.flow.settings.update(settings)
        self.flow.set_value_BC(pores=self.net.pores('back'), values=2)
        self.flow.set_value_BC(pores=self.net.pores('front'), values=0)

    def test_stokes_flow_power_law(self):
        hyd_cond = op.models.physics.hydraulic_conductance
        propname = 'throat.nonNewtonian_hydraulic_conductance'
        self.mod = hyd_cond.hagen_poiseuille_power_law
        self.phys.add_model(model=self.mod, regen_mode='normal',
                            propname=propname)
        self.phys.regenerate_models()
        self.flow.run()
        x = [0.,      0.,      0.,      0.,      0.,
             0.96051, 0.58898, 0.59932, 0.25297, 0.12773,
             1.30402, 1.06207, 1.3037,  1.05894, 0.44354,
             1.61413, 1.55325, 1.6475, 1.82346, 1.39743,
             2.,      2.,      2.,      2.,      2.]
        y = np.around(self.flow['pore.pressure'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = StokesFlowPowerLawTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
