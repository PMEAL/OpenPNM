import openpnm as op
import scipy as sp
from numpy.testing import assert_approx_equal
from numpy.testing import assert_allclose


class DiffusiveConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.diameter'] = 1
        self.geo['throat.diameter'] = 0.5
        self.geo['pore.area'] = 1
        self.geo['throat.area'] = 1
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.diffusivity'] = 1
        self.phase['pore.molecular_weight'] = 0.029
        self.phase['pore.temperature'] = 345
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)

    def test_ordinary_diffusion(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.15
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.25
        mod = op.models.physics.diffusive_conductance.ordinary_diffusion
        self.phys.add_model(propname='throat.o_diffusive_conductance',
                            model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.o_diffusive_conductance'].mean()
        assert_approx_equal(actual, desired=1.0)

    def test_ordinary_diffusion_with_zero_length_throats(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.15
        self.geo['throat.conduit_lengths.throat'] = 0.0
        self.geo['throat.conduit_lengths.pore2'] = 0.25
        mod = op.models.physics.diffusive_conductance.ordinary_diffusion
        self.phys.add_model(propname='throat.o_diffusive_conductance',
                            model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.o_diffusive_conductance'].mean()
        assert_approx_equal(actual, desired=2.5)

    def test_mixed_diffusion(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.15
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.25
        mod = op.models.physics.diffusive_conductance.mixed_diffusion
        self.phys.add_model(propname='throat.m_diffusive_conductance',
                            model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.m_diffusive_conductance'].mean()
        assert_approx_equal(actual, desired=0.9905265161584661)

    def test_taylor_aris_diffusion(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.15
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.25
        self.phase['pore.pressure'] = sp.linspace(0, 20, self.net.Np)
        self.phase['throat.hydraulic_conductance'] = 1
        mod = op.models.physics.diffusive_conductance.taylor_aris_diffusion
        self.phys.add_model(propname='throat.ta_diffusive_conductance',
                            model=mod)
        self.phys.regenerate_models()
        actual = sp.array([
            self.phys['throat.ta_diffusive_conductance'].mean(),
            self.phys['throat.ta_diffusive_conductance'].max(),
            self.phys['throat.ta_diffusive_conductance'].min()])
        desired = sp.array([1.03744, 1.10782, 1.00017])
        actual = sp.around(actual, decimals=5)
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
        assert sp.allclose(a=self.phys['throat.conductance'][0],
                           b=0.00084552)

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
        assert sp.allclose(a=self.phys['throat.conductance'][0],
                           b=0.00084552)


if __name__ == '__main__':

    t = DiffusiveConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
