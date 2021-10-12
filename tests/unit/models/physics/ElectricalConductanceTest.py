import numpy as np
import openpnm as op
from numpy.testing import assert_allclose
import openpnm.models.geometry.conduit_lengths as _conduit_lengths


class ElectricalConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[4, 4, 4])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.diameter'] = 1
        self.geo['throat.diameter'] = 1
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.electrical_conductivity'] = 1
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)

    def test_electrical_conductance_with_zero_length_throats(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.25
        self.geo['throat.conduit_lengths.throat'] = 0.0
        self.geo['throat.conduit_lengths.pore2'] = 0.15
        mpo = op.models.geometry.diffusive_size_factors.spheres_and_cylinders
        self.geo.add_model(propname='throat.diffusive_size_factors_zl',
                           model=mpo)
        mod = op.models.physics.electrical_conductance.generic_electrical
        self.phys.add_model(propname='throat.electrical_conductance', model=mod,
                            size_factors='throat.diffusive_size_factors_zl')
        self.phys.regenerate_models()
        actual = self.phys['throat.electrical_conductance'].mean()
        assert_allclose(actual, desired=0.785398, rtol=1e-5)

    def test_generic_electrical(self):
        self.geo['pore.diameter'] = 1.12
        self.geo['throat.diameter'] = 0.56
        L1, Lt, L2 = _conduit_lengths.spheres_and_cylinders(self.geo).T
        self.geo['throat.conduit_lengths.pore1'] = L1
        self.geo['throat.conduit_lengths.throat'] = Lt
        self.geo['throat.conduit_lengths.pore2'] = L2
        mpo2 = op.models.geometry.diffusive_size_factors.spheres_and_cylinders
        self.geo.add_model(propname='throat.diffusive_size_factors',
                           model=mpo2)
        mod2 = op.models.physics.electrical_conductance.generic_electrical
        self.phys.add_model(propname='throat.electrical_conductance_generic',
                            model=mod2)
        self.phys.regenerate_models()
        actual = np.mean(self.phys['throat.electrical_conductance_generic'])
        assert_allclose(actual, desired=0.61263, rtol=1e-5)


if __name__ == '__main__':

    t = ElectricalConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
