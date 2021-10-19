import openpnm as op
import numpy as _np
from numpy.testing import assert_allclose
import openpnm.models.geometry.conduit_lengths as _conduit_lengths


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
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.viscosity'] = 1e-5
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)

    def teardown_class(self):
        mgr = op.Workspace()
        mgr.clear()

    def test_generic_hydraulic(self):
        # Pass size factors as dict
        self.geo['throat.hydraulic_size_factors'] = {
            "pore1": 0.123, "throat": 0.981, "pore2": 0.551
        }
        mod = op.models.physics.hydraulic_conductance.generic_hydraulic
        self.phys.add_model(propname='throat.g_hydraulic_conductance', model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.g_hydraulic_conductance'].mean()
        assert_allclose(actual, desired=9120.483231751232)
        # Pass size factors as an array
        for elem in ["pore1", "throat", "pore2"]:
            del self.geo[f"throat.hydraulic_size_factors.{elem}"]
        self.geo['throat.hydraulic_size_factors'] = 0.896
        self.phys.regenerate_models("throat.g_hydraulic_conductance")
        actual = self.phys['throat.g_hydraulic_conductance'].mean()
        assert_allclose(actual, desired=89600.0)

    def test_hagen_poiseuille(self):
        del self.geo['throat.hydraulic_size_factors']
        self.geo['throat.hydraulic_size_factors'] = {
            "pore1": 0.123, "throat": 0.981, "pore2": 0.551
        }
        mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
        self.phys.add_model(propname='throat.hydraulic_conductance', model=mod)
        actual = self.phys['throat.hydraulic_conductance'].mean()
        assert_allclose(actual, desired=9120.483231751232)

    def test_hagen_poiseuille_2d(self):
        mod = op.models.geometry.hydraulic_size_factors.circles_and_rectangles
        self.geo.add_model(propname='throat.hydraulic_size_factor_2d', model=mod)
        mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
        self.phys.add_model(propname='throat.HP_hydraulic_conductance', model=mod,
                            size_factors='throat.hydraulic_size_factor_2d')
        actual = self.phys['throat.HP_hydraulic_conductance'].mean()
        assert_allclose(actual, desired=2972.1064015105035)

    def test_hagen_poiseuille_zero_length_throat(self):
        self.geo['throat.hydraulic_size_factors_zl'] = {
            "pore1": 0.123, "throat": _np.inf, "pore2": 0.551
        }
        mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
        self.phys.add_model(propname='throat.hydraulic_conductance',
                            model=mod,
                            size_factors='throat.hydraulic_size_factor_zl')
        actual = self.phys['throat.hydraulic_conductance'].mean()
        assert_allclose(actual, desired=9120.4832317)

    def test_valvatne_blunt(self):
        L1, Lt, L2 = _conduit_lengths.spheres_and_cylinders(self.geo).T
        self.geo['throat.conduit_lengths.pore1'] = L1
        self.geo['throat.conduit_lengths.throat'] = Lt
        self.geo['throat.conduit_lengths.pore2'] = L2
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.viscosity'] = 1e-5
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        mod = op.models.physics.hydraulic_conductance.valvatne_blunt
        sf = _np.sqrt(3) / 36.0
        self.geo['pore.shape_factor'] = _np.ones(self.geo.Np) * sf
        self.geo['throat.shape_factor'] = _np.ones(self.geo.Nt) * sf
        self.phys.add_model(propname='throat.valvatne_conductance', model=mod)
        actual = self.phys['throat.valvatne_conductance'].mean()
        desired = 2059.13571716
        assert_allclose(actual, desired=desired)


if __name__ == '__main__':

    t = HydraulicConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
