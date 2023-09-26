import numpy.testing as nptest

import openpnm as op


class WaterTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.phase = op.phase.Water(network=self.net)

    def test_generic_prop_with_temperature(self):
        T1 = 298.0
        T2 = 333.0
        vals = {
            'pore.viscosity':
                {
                    T1: 0.000893190947459112,
                    T2: 0.00046734380929754017,
                },
            'pore.density':
                {
                    T1: 996.9522269370573,
                    T2: 983.3169657125407,
                },
            'pore.thermal_conductivity':
                {
                    T1: 0.6104761069075901,
                    T2: 0.6500614873117142,
                },
            'pore.molar_density':
                {
                    T1: 55339.25794864455,
                    T2: 54582.3859364129,
                },
            'pore.surface_tension':
                {
                    T1: 0.07199532948823899,
                    T2: 0.06626423376953479,
                },
            'pore.vapor_pressure':
                {
                    T1: 3150.408314323942,
                    T2: 19812.906299267815,
                },
        }
        for prop in vals.keys():
            print(f'Testing {prop}')
            self.phase['pore.temperature'] = T1
            self.phase.regenerate_models()
            val = self.phase[prop][0]
            nptest.assert_allclose(val, vals[prop][T1], rtol=1e-10)
            self.phase['pore.temperature'] = T2
            self.phase.regenerate_models()
            val = self.phase[prop][0]
            nptest.assert_allclose(val, vals[prop][T2], rtol=1e-10)


if __name__ == '__main__':

    t = WaterTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()
