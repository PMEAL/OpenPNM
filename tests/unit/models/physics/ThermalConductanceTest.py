import openpnm as op
from numpy.testing import assert_approx_equal, assert_allclose
import openpnm.models.geometry.conduit_lengths as _conduit_lengths


class ThermalConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[4, 4, 4])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.geo['pore.diameter'] = 1.0
        self.geo['pore.area'] = 1.0
        self.geo['throat.diameter'] = 1.0
        self.geo['throat.length'] = 1e-9
        self.geo['throat.area'] = 1

    def test_thermal_conductance(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.25
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.15
        self.phase['pore.thermal_conductivity'] = 1
        mod = op.models.physics.thermal_conductance.series_resistors
        self.phys.add_model(propname='throat.thermal_conductance',
                            model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.thermal_conductance'].mean()
        assert_approx_equal(actual, desired=1.0)

    def test_thermal_conductance_with_zero_length_throats(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.25
        self.geo['throat.conduit_lengths.throat'] = 0.0
        self.geo['throat.conduit_lengths.pore2'] = 0.15
        self.phase['pore.thermal_conductivity'] = 1
        mod = op.models.physics.thermal_conductance.series_resistors
        self.phys.add_model(propname='throat.thermal_conductance',
                            model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.thermal_conductance'].mean()
        assert_approx_equal(actual, desired=2.5)
    
    def test_thermal_conductance_generic(self):
        self.geo['pore.diameter'] = 1.12
        self.geo['throat.diameter'] = 0.56
        L1, Lt, L2 = _conduit_lengths.spheres_and_cylinders(self.geo).T
        self.geo['throat.conduit_lengths.pore1'] = L1
        self.geo['throat.conduit_lengths.throat'] = Lt
        self.geo['throat.conduit_lengths.pore2'] = L2
        # old series resistors model, shape factor
        mpo = op.models.physics.poisson_shape_factors.ball_and_stick
        self.phys.add_model(propname="throat.poisson_shape_factors", model=mpo)
        mod1 = op.models.physics.thermal_conductance.series_resistors
        self.phys.add_model(propname='throat.thermal_conductance_from_mod',
                            model=mod1)
        self.phys.regenerate_models()
        # new series resistors model, size factor
        mpo2 = op.models.geometry.diffusive_size_factors.spheres_and_cylinders
        self.geo.add_model(propname="throat.diffusive_size_factors",
                           model=mpo2)
        mod2 = op.models.physics.thermal_conductance.series_resistors_generic
        self.phys.add_model(propname='throat.thermal_conductance_generic',
                            model=mod2)
        self.phys.regenerate_models()
        actual = self.phys['throat.thermal_conductance_from_mod']
        desired = self.phys['throat.thermal_conductance_generic']
        assert_allclose(actual, desired, rtol=1e-5)


if __name__ == '__main__':

    t = ThermalConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
