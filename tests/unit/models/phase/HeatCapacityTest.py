import openpnm as op
from numpy.testing import assert_allclose
import chemicals


class HeatCapacityTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])

    def test_generic_chemicals_for_pure_gas(self):
        mods = [
            # chemicals.heat_capacity.TRCCp,  # Requires constants
            # chemicals.heat_capacity.Shomate,  # Requires constants
            # chemicals.heat_capacity.Lastovka_Shaw,  # Requires similarity_valiable
            # chemicals.heat_capacity.Rowlinson_Poling,  # Requires Cpgm
            # chemicals.heat_capacity.Rowlinson_Bondi,  # Requires Cpgm
        ]
        n2 = op.phase.Species(network=self.net, species='nitrogen')
        vals = []
        for f in mods:
            vals.append(op.models.phase.chemicals_pure_prop(target=n2, f=f).mean())
        assert_allclose(vals, 0, rtol=.3)

    def test_generic_chemicals_for_pure_liq(self):
        mods = [
            # chemicals.heat_capacity.Zabransky_quasi_polynomial,  # Requires constants
            # chemicals.heat_capacity.Zabransky_cubic,  # Requires constants
            # chemicals.heat_capacity.Rowlinson_Poling,  # Requires constants
            # chemicals.heat_capacity.Rowlinson_Bondi,  # Requires Cpgm
            # chemicals.heat_capacity.Dadgostar_Shaw,  # Needs similarity_variable
            # chemicals.heat_capacity.Shomate,  # Requires constants
        ]
        h2o = op.phase.Species(network=self.net, species='water')
        vals = []
        for f in mods:
            vals.append(op.models.phase.chemicals_pure_prop(target=h2o, f=f).mean())
        assert_allclose(vals, 0, rtol=0.2)


if __name__ == '__main__':

    t = HeatCapacityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
