import openpnm as op
from numpy.testing import assert_approx_equal


class DiffusiveConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.diameter'] = 1.0
        self.geo['pore.area'] = 1.0
        self.geo['throat.diameter'] = 1.0
        self.geo['throat.length'] = 1e-9
        self.geo['throat.area'] = 1
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)

    def test_ordinary_diffusion(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.2
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.2
        self.geo['throat.equivalent_area.pore1'] = 0.2
        self.geo['throat.equivalent_area.throat'] = 0.2
        self.geo['throat.equivalent_area.pore2'] = 0.2
        mod = op.models.physics.diffusive_conductance.ordinary_diffusion
        self.phys.add_model(propname='throat.diffusive_conductance',
                            model=mod)
        self.phys.regenerate_models()
        assert_approx_equal(4.135096e-06,
                            self.phys['throat.diffusive_conductance'].mean())


if __name__ == '__main__':

    t = DiffusiveConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
