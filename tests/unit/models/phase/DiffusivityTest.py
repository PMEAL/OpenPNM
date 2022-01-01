import openpnm as op
from numpy.testing import assert_approx_equal


class DiffusivityTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.phase = op.phase.GenericPhase(network=self.net)
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.pressure'] = 101325  # Pa
        self.phase['pore.viscosity'] = 1.75e-5  # Pa.s

    def test_fuller(self):
        f = op.models.phase.diffusivity.fuller
        self.phase.add_model(propname='pore.diffusivity',
                             model=f,
                             MA=0.032,
                             MB=0.028,
                             vA=16.6,
                             vB=17.9)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.diffusivity'].mean(), 2.06754784e-05)

    def test_fuller_scaling(self):
        f = op.models.phase.diffusivity.fuller_scaling
        self.phase.add_model(propname='pore.diffusivity',
                             model=f,
                             DABo=1.79712526e-05,
                             Po=100000,
                             To=273)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.diffusivity'].mean(), 2.06754784e-05)

    def test_tyn_calus(self):
        f = op.models.phase.diffusivity.tyn_calus
        self.phase.add_model(propname='pore.diffusivity',
                             model=f,
                             VA=16.5,
                             VB=17.9,
                             sigma_A=1,
                             sigma_B=1)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.diffusivity'].mean(), 9.84851806e-05)

    def test_tyn_calus_scaling(self):
        f = op.models.phase.diffusivity.tyn_calus_scaling
        self.phase.add_model(propname='pore.diffusivity',
                             model=f,
                             DABo=5.26300839e-05,
                             mu_o=3e-5,
                             To=273)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.diffusivity'].mean(), 9.84851806e-05)


if __name__ == '__main__':

    t = DiffusivityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
