import openpnm as op
import numpy.testing as nptest


class MercuryTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.phase = op.phase.Mercury(network=self.net)

    def test_generic_prop_with_temperature(self):
        T1 = 298.0
        T2 = 333.0
        vals = {
            'pore.viscosity':
                {
                    T1: 0.0015360551647071996,
                    T2: 0.0014046786779692,
                },
            'pore.density':
                {
                    T1: 13544.82808,
                    T2: 13458.37668,
                },
            'pore.thermal_conductivity':
                {
                    T1: 8.514606495199999,
                    T2: 8.9719517682,
                },
            'pore.molar_density':
                {
                    T1: 67524.9418216262,
                    T2: 67093.95622912409,
                },
            'pore.surface_tension':
                {
                    T1: 0.4791000000000001,
                    T2: 0.46930000000000005,
                },
        }
        for prop in vals.keys():
            print(f'Testing {prop}')
            self.phase['pore.temperature'] = T1
            self.phase.regenerate_models()
            val = self.phase[prop][0]
            nptest.assert_array_almost_equal_nulp(val, vals[prop][T1])
            self.phase['pore.temperature'] = T2
            self.phase.regenerate_models()
            val = self.phase[prop][0]
            nptest.assert_array_almost_equal_nulp(val, vals[prop][T2])


if __name__ == '__main__':

    t = MercuryTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
