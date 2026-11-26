import openpnm as op
import numpy as np
from numpy.testing import assert_allclose
import openpnm.models.physics.hydraulic_conductance as hydraulic_conductance


class HydraulicConductanceTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.net['pore.diameter'] = 1.0
        self.net['throat.diameter'] = 0.5
        self.net['pore.area'] = 1.0
        self.net['throat.cross_sectional_area'] = 0.5
        self.phase = op.phase.Phase(network=self.net)
        self.phase['pore.viscosity'] = 1e-5
        self.size_factors = np.ones([self.net.Nt, 3])*(0.123, 0.981, 0.551)

    def teardown_class(self):
        mgr = op.Workspace()
        mgr.clear()

    def test_generic_hydraulic_size_factors_as_dict(self):
        self.net['throat.hydraulic_size_factors'] = self.size_factors
        mod = op.models.physics.hydraulic_conductance.generic_hydraulic
        self.phase.add_model(propname='throat.g_hydraulic_conductance', model=mod)
        self.phase.regenerate_models()
        actual = self.phase['throat.g_hydraulic_conductance'].mean()
        assert_allclose(actual, desired=9120.483231751232)
        del self.net["throat.hydraulic_size_factors"]
        # del self.phase['throat.g_hydraulic_conductance']

    def test_generic_hydraulic_size_factors_as_array(self):
        self.net['throat.hydraulic_size_factors'] = 0.896
        self.phase.regenerate_models()
        actual = self.phase['throat.g_hydraulic_conductance'].mean()
        assert_allclose(actual, desired=89600.0)
        del self.net["throat.hydraulic_size_factors"]
        # del self.phase['throat.g_hydraulic_conductance']

    def test_hagen_poiseuille(self):
        self.net['throat.hydraulic_size_factors'] = self.size_factors
        mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
        self.phase.add_model(propname='throat.hydraulic_conductance', model=mod)
        actual = self.phase['throat.hydraulic_conductance'].mean()
        assert_allclose(actual, desired=9120.483231751232)

    def test_valvatne_blunt(self):
        self.net['throat.conduit_lengths'] = np.ones([self.net.Nt, 3])*(0.15, 0.2, 0.26)
        self.net['pore.shape_factor'] = np.ones(self.net.Np) * 0.15
        self.net['throat.shape_factor'] = np.ones(self.net.Nt) * 0.15
        mod = op.models.physics.hydraulic_conductance.valvatne_blunt
        self.phase.add_model(propname='throat.valvatne_conductance', model=mod)
        actual = self.phase['throat.valvatne_conductance'].mean()
        desired = 6198.347107
        assert_allclose(actual, desired=desired)
        
    def test_conductance_shape(self):
        available_models = [
            getattr(hydraulic_conductance, model_name)
            for model_name in dir(hydraulic_conductance)
            if callable(getattr(hydraulic_conductance, model_name))
            ]
        for model in available_models:
            G = model(phase=self.phase)
            assert_allclose(G.shape, (self.net.Nt, 2), rtol=0)

if __name__ == '__main__':

    t = HydraulicConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
