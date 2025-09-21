import numpy.testing as nptest

import openpnm as op


class AirTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.phase = op.phase.Air(network=self.net)

    def test_generic_prop_with_temperature(self):
        T1 = 298.0
        T2 = 333.0
        vals = {
            'pore.viscosity':
                {
                    T1: 1.8444445170912796e-05,
                    T2: 2.00725794239133e-05,
                },
            'pore.density':
                {
                    T1: 1.1798234085231065,
                    T2: 1.0558179451648222,
                },
            'pore.thermal_conductivity':
                {
                    T1: 0.026369425206800003,
                    T2: 0.028787674351300006,
                },
            'pore.diffusivity':
                {
                    T1: 2.0946939351898755e-05,
                    T2: 2.5440141454852606e-05,
                },
            'pore.molar_density':
                {
                    T1: 40.89461870781484,
                    T2: 36.596385510296765,
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

    t = AirTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()
