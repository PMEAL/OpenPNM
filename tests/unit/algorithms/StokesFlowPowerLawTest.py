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

    def test_stokes_flow_power_law_cylinder(self):
        hyd_cond = op.models.physics.hydraulic_conductance
        propname = 'throat.nonNewtonian_hydraulic_conductance'
        self.mod = hyd_cond.hagen_poiseuille_power_law
        self.phys.add_model(model=self.mod, regen_mode='normal',
                            propname=propname, shape='cylinder')
        self.phys.regenerate_models()
        self.flow.run()
        x = [0.,      0.,      0.,      0.,      0.,
             0.96144, 0.59069, 0.59915, 0.25231, 0.12742,
             1.30443, 1.06373, 1.30464, 1.05821, 0.44284,
             1.61415, 1.55411, 1.64821, 1.82377, 1.39799,
             2.,      2.,      2.,      2.,      2.]
        y = np.around(self.flow['pore.pressure'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_stokes_flow_power_law_cone(self):
        hyd_cond = op.models.physics.hydraulic_conductance
        propname = 'throat.nonNewtonian_hydraulic_conductance'
        self.mod = hyd_cond.hagen_poiseuille_power_law
        self.phys.add_model(model=self.mod, regen_mode='normal',
                            propname=propname, shape='cone')
        self.phys.regenerate_models()
        self.flow.run()
        x = [0.,      0.,      0.,      0.,      0.,
             0.95375, 0.58708, 0.59716, 0.25443, 0.13113,
             1.29703, 1.05664, 1.29892, 1.05618, 0.44736,
             1.61168, 1.5502,  1.6427, 1.81941, 1.39755,
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
