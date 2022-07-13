import openpnm as op
import chemicals
from thermo import Chemical
from numpy.testing import assert_approx_equal, assert_allclose


class ThermalConductivityTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])

    def test_water(self):
        f = op.models.phase.thermal_conductivity.water_correlation
        phase = op.phase.Species(network=self.net, species='water')
        phase.add_model(propname='pore.thermal_conductivity', model=f)
        assert_approx_equal(phase['pore.thermal_conductivity'].mean(),
                            0.61047611)

    def test_chung(self):
        f = op.models.phase.thermal_conductivity.gas_pure_chung
        phase = op.phase.Species(network=self.net, species='nitrogen')
        a = Chemical('nitrogen')
        phase['pore.heat_capacity'] = a.Cpg
        phase.add_model(propname='pore.viscosity',
                        model=op.models.phase.chemicals_pure_prop_wrapper,
                        f=chemicals.viscosity.viscosity_gas_Gharagheizi)
        phase.add_model(propname='pore.thermal_conductivity', model=f)
        phase.regenerate_models()
        assert_approx_equal(phase['pore.thermal_conductivity'].mean(),
                            0.82491818)

    def test_sato(self):
        f = op.models.phase.thermal_conductivity.liquid_pure_sato_riedel
        phase = op.phase.Species(network=self.net, species='water')
        phase.add_model(propname='pore.thermal_conductivity', model=f)
        phase.regenerate_models()
        assert_approx_equal(phase['pore.thermal_conductivity'].mean(),
                            0.297730167)

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
            vals.append(op.models.phase.chemicals_wrapper(target=a, f=f).mean())
        assert_allclose(vals, 2.898e-2, rtol=.2)

    def test_generic_chemicals_for_pure_liq(self):
        mods = [
            chemicals.thermal_conductivity.Bahadori_liquid,
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
            vals.append(op.models.phase.chemicals_wrapper(target=h2o, f=f).mean())
        assert_allclose(vals, 2.898e-1, rtol=1.5)


if __name__ == '__main__':

    t = ThermalConductivityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
