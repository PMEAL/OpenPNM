import openpnm as op
import numpy as _np
from numpy.testing import assert_allclose
from openpnm.utils import remove_prop_deep


class HydraulicConductanceTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.diameter'] = 1.0
        self.geo['throat.diameter'] = 0.5
        self.geo['pore.area'] = 1.0
        self.geo['throat.area'] = 0.5
        self.phase = op.phase.GenericPhase(network=self.net)
        self.phase['pore.viscosity'] = 1e-5
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.size_factors_dict = {"pore1": 0.123, "throat": 0.981, "pore2": 0.551}

    def teardown_class(self):
        mgr = op.Workspace()
        mgr.clear()

    def test_generic_hydraulic_size_factors_as_dict(self):
        self.geo['throat.hydraulic_size_factors'] = self.size_factors_dict
        mod = op.models.physics.hydraulic_conductance.generic_hydraulic
        self.phys.add_model(propname='throat.g_hydraulic_conductance', model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.g_hydraulic_conductance'].mean()
        assert_allclose(actual, desired=9120.483231751232)
        remove_prop_deep(self.geo, "throat.hydraulic_size_factors")

    def test_generic_hydraulic_size_factors_as_array(self):
        self.geo['throat.hydraulic_size_factors'] = 0.896
        self.phys.regenerate_models("throat.g_hydraulic_conductance")
        actual = self.phys['throat.g_hydraulic_conductance'].mean()
        assert_allclose(actual, desired=89600.0)
        remove_prop_deep(self.geo, "throat.hydraulic_size_factors")

    def test_hagen_poiseuille(self):
        self.geo['throat.hydraulic_size_factors'] = self.size_factors_dict
        mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
        self.phys.add_model(propname='throat.hydraulic_conductance', model=mod)
        actual = self.phys['throat.hydraulic_conductance'].mean()
        assert_allclose(actual, desired=9120.483231751232)

    def test_valvatne_blunt(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.15
        self.geo['throat.conduit_lengths.throat'] = 0.2
        self.geo['throat.conduit_lengths.pore2'] = 0.26
        self.geo['pore.shape_factor'] = _np.ones(self.geo.Np) * 0.15
        self.geo['throat.shape_factor'] = _np.ones(self.geo.Nt) * 0.15
        mod = op.models.physics.hydraulic_conductance.valvatne_blunt
        self.phys.add_model(propname='throat.valvatne_conductance', model=mod)
        actual = self.phys['throat.valvatne_conductance'].mean()
        desired = 6198.347107
        assert_allclose(actual, desired=desired)


if __name__ == '__main__':

    t = HydraulicConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
