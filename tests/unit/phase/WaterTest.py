import openpnm as op
import numpy.testing as nptest


class WaterTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.water = op.phase.Water(network=self.net)

    def test_change_temperature(self):
        assert self.water['pore.temperature'][0] == 298.0
        self.water.regenerate_models()
        mu = self.water['pore.viscosity'][0]
        nptest.assert_array_almost_equal_nulp(mu, 0.000893190947459112)
        self.water['pore.temperature'][0] = 333.0
        assert self.water['pore.temperature'][0] == 333.0
        self.water.regenerate_models()
        mu = self.water['pore.viscosity'][0]
        nptest.assert_array_almost_equal_nulp(mu, 0.00046734380929754017)


if __name__ == '__main__':

    t = WaterTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
