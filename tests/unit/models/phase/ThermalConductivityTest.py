import openpnm as op
from numpy.testing import assert_approx_equal, assert_array_almost_equal, assert_allclose
import chemicals


class ThermalConductivityTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.phase = op.phase.Phase(network=self.net)
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.salinity'] = 0.0  # g/kg
        self.phase['pore.viscosity'] = 0.001  # Pa.s
        self.phase['pore.critical_temperature'] = 647.15  # K
        self.phase['pore.molecular_weight'] = 0.018  # kg/mol
        self.phase['pore.boiling_point'] = 373.15  # K
        self.phase['pore.heat_capacity'] = 75.28  # J/mol K
        self.phase['pore.acentric_factor'] = 11.5  # J/mol K

    def test_water(self):
        f = op.models.phase.thermal_conductivity.water
        self.phase.add_model(propname='pore.thermal_conductivity',
                             model=f)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.thermal_conductivity'].mean(),
                            0.61047611)

    def test_chung(self):
        f = op.models.phase.thermal_conductivity.chung
        self.phase.add_model(propname='pore.thermal_conductivity',
                             model=f)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.thermal_conductivity'].mean(),
                            0.62063913)

    def test_sato(self):
        f = op.models.phase.thermal_conductivity.sato
        self.phase.add_model(propname='pore.thermal_conductivity',
                             model=f)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.thermal_conductivity'].mean(),
                            0.29787023)

    def test_generic_chemicals_for_pure_gas(self):
        mods = [
            chemicals.thermal_conductivity.Bahadori_gas,
            chemicals.thermal_conductivity.Gharagheizi_gas,
            # chemicals.thermal_conductivity.Eli_Hanley,  # Needs Cvm
            # chemicals.thermal_conductivity.Chung,  # Needs Cvm
            # chemicals.thermal_conductivity.DIPPR9B,  # Needs Cvm
            # chemicals.thermal_conductivity.Eucken_modified,  # Needs Cvm
            # chemicals.thermal_conductivity.Eucken,  # Needs Cvm
            # chemicals.thermal_conductivity.Stiel_Thodos_dense,  # Needs Vm
            # chemicals.thermal_conductivity.Eli_Hanley_dense,  # Needs Cvm
            # chemicals.thermal_conductivity.Chung_dense,  # Needs Cvm
        ]
        a = op.phase.Species(network=self.net, species='nitrogen')
        vals = []
        for f in mods:
            print(f)
            vals.append(op.models.phase.chemicals_pure_prop(target=a, f=f).mean())
        assert_allclose(vals, 2.898e-1, rtol=2)

    def test_generic_chemicals_for_pure_liq(self):
        mods = [
            chemicals.thermal_conductivity.Sheffy_Johnson,
            chemicals.thermal_conductivity.Sato_Riedel,
            chemicals.thermal_conductivity.Lakshmi_Prasad,
            chemicals.thermal_conductivity.Gharagheizi_liquid,
            # chemicals.thermal_conductivity.Nicola_original,  # Needs Hfus
            chemicals.thermal_conductivity.Nicola,
            chemicals.thermal_conductivity.Bahadori_liquid,
            # chemicals.thermal_conductivity.DIPPR9G,  # Needs kl
            # chemicals.thermal_conductivity.Missenard,  # Needs kl
        ]
        h2o = op.phase.Species(network=self.net, species='water')
        vals = []
        for f in mods:
            vals.append(op.models.phase.chemicals_pure_prop(target=h2o, f=f).mean())
        assert_allclose(vals, 2.898e-1, rtol=2)


if __name__ == '__main__':

    t = ThermalConductivityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
