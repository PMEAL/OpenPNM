import numpy as np
import openpnm as op
from openpnm.phases import mixtures
from numpy.testing import assert_allclose


class IonicConductanceTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.diameter'] = 1.12
        self.geo['throat.diameter'] = 0.56
        self.geo['pore.area'] = 1.0
        self.geo['throat.cross_sectional_area'] = 1.0
        self.geo['pore.volume'] = 1.0
        self.geo['throat.volume'] = 1.0
        self.geo["throat.diffusive_size_factors"] = {
            "pore1": 0.123, "throat": 0.981, "pore2": 0.551
        }
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.permittivity'] = 78.0
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)

    def test_generic_ionic_poisson(self):
        mod = op.models.physics.ionic_conductance.poisson
        self.phys.add_model(propname='throat.ionic_conductance', model=mod)
        self.phys.regenerate_models()
        actual = np.mean(self.phys['throat.ionic_conductance'])
        assert_allclose(actual, desired=6.298849e-11, rtol=1e-5)

    def test_ionic_conductance_electroneutrality(self):
        sw = mixtures.SalineWater(network=self.net)
        Na = sw.components[f'Na_{sw.name}']
        Cl = sw.components[f'Cl_{sw.name}']
        self.phys = op.physics.GenericPhysics(network=self.net, phase=sw, geometry=self.geo)
        sw[f'pore.concentration.{Cl.name}'] = np.linspace(10, 20, 27) * 0.5
        sw[f'pore.concentration.{Na.name}'] = np.linspace(10, 20, 27)
        electroneutrality = op.models.physics.ionic_conductance.electroneutrality
        self.phys.add_model(propname='throat.ionic_conductance',
                            model=electroneutrality, ions=[Na.name, Cl.name])
        self.phys.regenerate_models()
        actual = np.mean(self.phys['throat.ionic_conductance'])
        assert_allclose(actual, desired=0.0116345, rtol=1e-5)


if __name__ == '__main__':

    t = IonicConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('Running test: '+item)
            t.__getattribute__(item)()
