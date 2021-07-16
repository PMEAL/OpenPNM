import openpnm as op
from numpy.testing import assert_allclose


class IonicConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[4, 4])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.diameter'] = 1.12
        self.geo['throat.diameter'] = 0.56
        self.geo['pore.area'] = 1
        self.geo['throat.area'] = 1
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.permittivity'] = 78
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)

    def test_ionic_conductance_poisson_2D(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.15
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.25
        mod = op.models.physics.ionic_conductance.poisson
        self.phys.add_model(propname='throat.ionic_conductance_2D_from_mod',
                            model=mod)
        self.phys.regenerate_models()
        self.phys.add_model(propname='throat.ionic_conductance_2D_input',
                            model=mod,
                            pore_area='pore.diameter',
                            throat_area='throat.diameter')
        self.phys.regenerate_models()
        actual = self.phys['throat.ionic_conductance_2D_from_mod']
        desired = self.phys['throat.ionic_conductance_2D_input']
        assert_allclose(actual, desired, rtol=1e-5)


if __name__ == '__main__':

    t = IonicConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
