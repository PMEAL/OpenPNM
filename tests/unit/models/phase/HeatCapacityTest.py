import chemicals
from numpy.testing import assert_allclose
from thermo import Chemical

import openpnm as op


class HeatCapacityTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])

    def test_generic_chemicals_for_pure_gas(self):
        mods = [
            # chemicals.heat_capacity.TRCCp,  # Requires constants
            # chemicals.heat_capacity.Shomate,  # Requires constants
            # chemicals.heat_capacity.Lastovka_Shaw,  # Requires similarity_valiable
        ]
        n2 = op.phase.Species(network=self.net, species='nitrogen')
        vals = []
        for f in mods:
            vals.append(op.models.phase.chemicals_wrapper(phase=n2, f=f).mean())
        assert_allclose(vals, 0, rtol=.3)

    def test_generic_chemicals_for_pure_liq(self):
        mods = [
            # chemicals.heat_capacity.Zabransky_quasi_polynomial,  # Requires constants
            # chemicals.heat_capacity.Zabransky_cubic,  # Requires constants
            # chemicals.heat_capacity.Rowlinson_Poling,  # Requires Cpgm
            # chemicals.heat_capacity.Rowlinson_Bondi,  # Requires Cpgm
            # chemicals.heat_capacity.Dadgostar_Shaw,  # Needs similarity_variable
            # chemicals.heat_capacity.Shomate,  # Requires constants
        ]
        h2o = op.phase.Species(network=self.net, species='water')
        a = Chemical('h2o')
        h2o['pore.heat_capacity_gas'] = a.Cpgm
        vals = []
        for f in mods:
            vals.append(op.models.phase.chemicals_wrapper(phase=h2o, f=f).mean())
        assert_allclose(vals, 0, rtol=0.2)

    def test_custom_implementations(self):
        pn = op.network.Demo()
        ch4 = op.phase.Species(network=pn, species='ch4')
        ch4['pore.temperature'] = 100
        a = Chemical('ch4')
        a.T = 100
        ch4['pore.heat_capacity_gas'] = a.Cpgm
        Cp_ref = chemicals.heat_capacity.Rowlinson_Poling(
            T=ch4['pore.temperature'][0],
            Tc=ch4['param.critical_temperature'],
            omega=ch4['param.acentric_factor'],
            Cpgm=ch4['pore.heat_capacity_gas'][0])
        Cp_calc = op.models.phase.heat_capacity.liquid_pure_rp(ch4)
        assert_allclose(Cp_ref, Cp_calc, rtol=1e-10)


if __name__ == '__main__':

    t = HeatCapacityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()
