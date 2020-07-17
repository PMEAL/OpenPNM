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
        self.phase['pore.viscosity'] = 0.0011

        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.flow = op.algorithms.StokesFlow(network=self.net,
                                             phase=self.phase)
        settings = {'cache_A': False, 'cache_b': False,
                    'solver_type': 'spsolve', 'max_iter': 200,
                    'solver_tol': 1e-6}
        self.flow.settings.update(settings)
        self.flow.set_value_BC(pores=self.net.pores('back'), values=2)
        self.flow.set_value_BC(pores=self.net.pores('front'), values=0)

    def test_stokes_flow_cylinder(self):
        hyd_cond = op.models.physics.hydraulic_conductance
        propname = 'throat.hydraulic_conductance'
        self.mod = hyd_cond.hagen_poiseuille
        self.phys.add_model(model=self.mod, regen_mode='normal',
                            propname=propname, shape='cylinder')
        self.phys.regenerate_models()
        self.flow.run()
        x = [0.,      0.,      0.,      0.,      0.,
             1.04615, 0.73872, 0.57359, 0.16617, 0.08652,
             1.34667, 1.20966, 1.37768, 1.00439, 0.35137,
             1.62217, 1.64073, 1.71774, 1.85626, 1.48442,
             2.,      2.,      2.,      2.,      2.]
        y = np.around(self.flow['pore.pressure'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_stokes_flow_cone(self):
        hyd_cond = op.models.physics.hydraulic_conductance
        propname = 'throat.hydraulic_conductance'
        self.mod = hyd_cond.hagen_poiseuille
        self.phys.add_model(model=self.mod, regen_mode='normal',
                            propname=propname, shape='cone')
        self.phys.regenerate_models()
        self.flow.run()
        x = [0.,      0.,      0.,      0.,      0.,
             1.03009, 0.73073, 0.56994, 0.17297, 0.09413,
             1.33267, 1.19722, 1.36541, 0.99959, 0.36484,
             1.61667, 1.6334,  1.70666, 1.84742, 1.48111,
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
