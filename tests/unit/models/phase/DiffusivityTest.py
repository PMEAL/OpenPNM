import openpnm as op
from numpy.testing import assert_approx_equal


class DiffusivityTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.phase = op.phase.Phase(network=self.net)
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.pressure'] = 101325  # Pa
        self.phase['pore.viscosity'] = 1.75e-5  # Pa.s

    def test_fuller(self):
        f = op.models.phase.diffusivity.gas_mixture_fesg
        self.phase.add_model(propname='pore.diffusivity',
                             model=f,
                             MWs=[31.9988, 28.0134],
                             Vdms=[16.6, 17.9],)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.diffusivity'][0], 2.094693935e-05)


if __name__ == '__main__':

    t = DiffusivityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
