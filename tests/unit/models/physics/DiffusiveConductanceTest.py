import numpy as np
import openpnm as op
from numpy.testing import assert_allclose


class DiffusiveConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.diameter'] = 1.
        self.geo['throat.diameter'] = 0.5
        self.geo['pore.area'] = 1.
        self.geo['throat.area'] = 1.
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.diffusivity'] = 1.3
        self.phase['pore.molecular_weight'] = 0.029
        self.phase['pore.temperature'] = 345
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)

    def test_generic_diffusive(self):
        # Pass size factors as dict
        self.geo['throat.diffusive_size_factors'] = {
            "pore1": 0.123, "throat": 0.981, "pore2": 0.551
        }
        mod = op.models.physics.diffusive_conductance.generic_diffusive
        self.phys.add_model(propname='throat.g_diffusive_conductance', model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.g_diffusive_conductance'].mean()
        assert_allclose(actual, desired=0.091204832 * 1.3)
        # Pass size factors as an array
        for elem in ["pore1", "throat", "pore2"]:
            del self.geo[f"throat.diffusive_size_factors.{elem}"]
        self.geo['throat.diffusive_size_factors'] = 0.896
        self.phys.regenerate_models("throat.g_diffusive_conductance")
        actual = self.phys['throat.g_diffusive_conductance'].mean()
        assert_allclose(actual, desired=0.896 * 1.3)

    def test_generic_diffusive_partial_domain(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        geo = op.geometry.GenericGeometry(network=net,
                                          pores=net.Ps,
                                          throats=net.Ts[0:5])
        phase = op.phases.GenericPhase(network=net)
        phase['pore.diffusivity'] = 1.3
        geo['throat.diffusive_size_factors'] = 0.5
        phys = op.physics.GenericPhysics(network=net,
                                         phase=phase,
                                         geometry=geo)
        mod = op.models.physics.diffusive_conductance.generic_diffusive
        phys.add_model(propname='throat.g_diffusive_conductance', model=mod)
        phys.regenerate_models()
        actual = phys['throat.g_diffusive_conductance'].mean()
        assert (actual == 0.65 and len(phys['throat.g_diffusive_conductance']) == 5)

    def test_ordinary_diffusion(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.15
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.25
        mod = op.models.physics.diffusive_conductance.ordinary_diffusion
        self.phys.add_model(propname='throat.o_diffusive_conductance', model=mod)
        actual = self.phys['throat.o_diffusive_conductance'].mean()
        assert_allclose(actual, desired=1.3)

    def test_ordinary_diffusion_2d(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.15
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.25
        mod = op.models.physics.diffusive_conductance.ordinary_diffusion_2d
        self.phys.add_model(propname='throat.o_diffusive_conductance', model=mod)
        actual = self.phys['throat.o_diffusive_conductance'].mean()
        assert_allclose(actual, desired=0.8125)

    def test_ordinary_diffusion_with_zero_length_throats(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.15
        self.geo['throat.conduit_lengths.throat'] = 0.0
        self.geo['throat.conduit_lengths.pore2'] = 0.25
        mod = op.models.physics.diffusive_conductance.ordinary_diffusion
        self.phys.add_model(propname='throat.o_diffusive_conductance', model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.o_diffusive_conductance'].mean()
        assert_allclose(actual, desired=2.5 * 1.3)

    def test_mixed_diffusion(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.15
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.25
        mod = op.models.physics.diffusive_conductance.mixed_diffusion
        self.phys.add_model(propname='throat.m_diffusive_conductance', model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.m_diffusive_conductance'].mean()
        assert_allclose(actual, desired=1.284035187705028)

    def test_multiphase_diffusion(self):
        np.random.seed(50)
        net = op.network.Cubic(shape=[1, 6, 1])
        geom = op.geometry.SpheresAndCylinders(network=net,
                                               pores=net.Ps, throats=net.Ts)
        air = op.phases.Air(network=net)
        water = op.phases.Water(network=net)
        m = op.phases.MultiPhase(phases=[air, water], project=net.project)
        m.set_occupancy(phase=air, pores=[0, 1, 2])
        m.set_occupancy(phase=water, pores=[3, 4, 5])
        const = op.models.misc.constant
        K_water_air = 0.5
        m.set_binary_partition_coef(
            phases=[water, air], model=const, value=K_water_air
        )
        m._set_automatic_throat_occupancy()
        phys = op.physics.GenericPhysics(network=net, phase=m, geometry=geom)
        mdiff = op.models.physics.diffusive_conductance.multiphase_diffusion
        phys.add_model(propname="throat.diffusive_conductance", model=mdiff)
        g = phys["throat.diffusive_conductance"]
        # Diffusive conductance for MultiPhase must be Nt by 2 (not Nt by 1)
        assert g.shape == (net.Nt, 2)
        # Columns 1, 2 of conductance must be equal except for interface throat
        assert_allclose(g[:, 0][[0, 1, 3, 4]], g[:, 1][[0, 1, 3, 4]])
        # G12 and G21 at interface must differ (ratio must be K_water_air)
        assert_allclose(g[2, 0] / g[2, 1], 1 / K_water_air)
        assert_allclose(g.mean(), 2.139269316e-7)

    def test_taylor_aris_diffusion(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.15
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.25
        self.phase['pore.pressure'] = np.linspace(0, 20, self.net.Np)
        self.phase['throat.hydraulic_conductance'] = 1
        mod = op.models.physics.diffusive_conductance.taylor_aris_diffusion
        self.phys.add_model(propname='throat.ta_diffusive_conductance',
                            model=mod)
        self.phys.regenerate_models()
        actual = np.array([
            self.phys['throat.ta_diffusive_conductance'].mean(),
            self.phys['throat.ta_diffusive_conductance'].max(),
            self.phys['throat.ta_diffusive_conductance'].min()])
        desired = np.array([1.32879664, 1.38293965, 1.3001327])
        assert_allclose(actual, desired)

    def test_classic_ordinary_diffusion(self):
        self.geo['pore.diameter'] = 1.0
        self.geo['throat.diameter'] = 1.0
        self.geo['throat.length'] = 1e-9
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)
        mod = op.models.physics.diffusive_conductance.classic_ordinary_diffusion
        self.phys.add_model(propname='throat.conductance', model=mod)
        assert np.allclose(a=self.phys['throat.conductance'][0], b=0.00084552)

    def test_classic_ordinary_diffusion_with_zero_length_throats(self):
        self.geo['pore.diameter'] = 1.0
        self.geo['throat.diameter'] = 1.0
        self.geo['throat.length'] = 0.0
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)
        mod = op.models.physics.diffusive_conductance.classic_ordinary_diffusion
        self.phys.add_model(propname='throat.conductance', model=mod)
        assert np.allclose(a=self.phys['throat.conductance'][0], b=0.00084552)


if __name__ == '__main__':

    t = DiffusiveConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
