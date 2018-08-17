import openpnm as op
from numpy.testing import assert_approx_equal


class DiffusiveConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.area'] = 1
        self.geo['throat.area'] = 1
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.diffusivity'] = 1
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)

    def test_ordinary_diffusion(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.15
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.25
        mod = op.models.physics.diffusive_conductance.ordinary_diffusion
        self.phys.add_model(propname='throat.diffusive_conductance',
                            pore_diffusivity='pore.diffusivity',
                            model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.diffusive_conductance'].mean()
        assert_approx_equal(actual, desired=1.0)


if __name__ == '__main__':

    t = DiffusiveConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
