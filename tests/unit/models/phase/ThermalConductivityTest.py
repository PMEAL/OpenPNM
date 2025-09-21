import inspect

import chemicals
import numpy as np
from numpy.testing import assert_allclose, assert_approx_equal
from thermo import Chemical

import openpnm as op
from openpnm.utils import get_mixture_model_args


class ThermalConductivityTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])

    def test_water_correlation(self):
        f = op.models.phase.thermal_conductivity.water_correlation
        phase = op.phase.Species(network=self.net, species='water')
        phase.add_model(propname='pore.thermal_conductivity', model=f)
        assert_approx_equal(phase['pore.thermal_conductivity'].mean(),
                            0.61047611)

    def test_pure_gas(self):
        pn = op.network.Demo()
        gas = op.phase.Species(network=pn, species='ch4')
        gas.add_model(propname='pore.thermal_conductivity',
                      model=op.models.phase.thermal_conductivity.gas_pure_gismr)
        temp = chemicals.thermal_conductivity.Gharagheizi_gas(
            T=gas['pore.temperature'][0],
            Tb=gas['param.boiling_temperature'],
            Pc=gas['param.critical_pressure'],
            MW=gas['param.molecular_weight'],
            omega=gas['param.acentric_factor'],
        )
        assert_allclose(temp, gas['pore.thermal_conductivity'][0], rtol=1e-10)

    def test_pure_liquid(self):
        pn = op.network.Demo()
        liq = op.phase.Species(network=pn, species='h2o')
        liq.add_model(propname='pore.thermal_conductivity',
                      model=op.models.phase.thermal_conductivity.liquid_pure_gismr)
        temp = chemicals.thermal_conductivity.Gharagheizi_liquid(
            T=liq['pore.temperature'][0],
            Tb=liq['param.boiling_temperature'],
            Pc=liq['param.critical_pressure'],
            MW=liq['param.molecular_weight'],
            omega=liq['param.acentric_factor'],
        )
        assert_allclose(temp, liq['pore.thermal_conductivity'][0], rtol=1e-10)

    def test_gas_mixture(self):
        pn = op.network.Demo()
        gas1 = op.phase.Species(network=pn, species='ch4')
        gas2 = op.phase.Species(network=pn, species='co2')
        gas1.add_model(propname='pore.thermal_conductivity',
                       model=op.models.phase.thermal_conductivity.gas_pure_gismr)
        gas2.add_model(propname='pore.thermal_conductivity',
                       model=op.models.phase.thermal_conductivity.gas_pure_gismr)
        mix = op.phase.GasMixture(network=pn, components=[gas1, gas2])
        mix.y(gas1, 0.5)
        mix.y(gas2, 0.5)
        mix.add_model(propname='pore.thermal_conductivity',
                      model=op.models.phase.thermal_conductivity.gas_mixture_whz)
        args = get_mixture_model_args(mix, composition='zs',
                                      args={'ks': 'pore.thermal_conductivity',
                                            'MWs': 'param.molecular_weight'})
        temp = chemicals.thermal_conductivity.Wassiljewa_Herning_Zipperer(**args)
        assert_allclose(temp, mix['pore.thermal_conductivity'][0], rtol=1e-10)

    def test_liquid_mixture(self):
        pn = op.network.Demo()
        liq1 = op.phase.Species(network=pn, species='h2o')
        liq2 = op.phase.Species(network=pn, species='etoh')
        liq1.add_model(propname='pore.thermal_conductivity',
                       model=op.models.phase.thermal_conductivity.liquid_pure_gismr)
        liq2.add_model(propname='pore.thermal_conductivity',
                       model=op.models.phase.thermal_conductivity.liquid_pure_gismr)
        mix = op.phase.LiquidMixture(network=pn, components=[liq1, liq2])
        mix.x(liq1, 0.5)
        mix.x(liq2, 0.5)
        mix.add_model(propname='pore.thermal_conductivity',
                      model=op.models.phase.thermal_conductivity.liquid_mixture_DIPPR9H)
        args = get_mixture_model_args(mix, composition='ws',
                                      args={'ks': 'pore.thermal_conductivity'})
        temp = chemicals.thermal_conductivity.DIPPR9H(**args)
        assert_allclose(temp, mix['pore.thermal_conductivity'][0], rtol=1e-10)

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
            vals.append(op.models.phase.chemicals_wrapper(phase=a, f=f).mean())
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
            vals.append(op.models.phase.chemicals_wrapper(phase=h2o, f=f).mean())
        assert_allclose(vals, 2.898e-1, rtol=1.5)

    def test_pure_liquid_thermal_conductivity_models(self):
        pn = op.network.Demo()
        phase = op.phase.Species(network=pn, species='water')
        argmap = op.models.phase.default_argmap
        vals = {
            'Gharagheizi_liquid': op.models.phase.thermal_conductivity.liquid_pure_gismr,
            'Sato_Riedel': op.models.phase.thermal_conductivity.liquid_pure_sr
        }
        for k, v in vals.items():
            print(f'Testing {k}')
            f = getattr(chemicals.thermal_conductivity, k)
            args = inspect.getfullargspec(f)[0]
            kwargs = {i: np.atleast_1d(phase[argmap[i]])[0] for i in args}
            ref = f(**kwargs)
            val = v(phase)[0]
            assert_allclose(ref, val, rtol=1e-13)

    def test_liquid_mixture_thermal_conductivity_models(self):
        pn = op.network.Demo()
        a = op.phase.StandardLiquid(network=pn, species='water')
        b = op.phase.StandardLiquid(network=pn, species='ethanol')
        phase = op.phase.StandardLiquidMixture(network=pn, components=[a, b])
        phase.add_model(propname='pore.mass_fraction',
                        model=op.models.phase.mixtures.mole_to_mass_fraction,
                        regen_mode='deferred')
        phase.x(a, 0.5)
        phase.x(b, 0.5)
        phase.regenerate_models()
        argmap = op.models.phase.default_argmap
        argmap.update(op.models.phase.mixture_argmap)
        vals = {
            'DIPPR9H': op.models.phase.thermal_conductivity.liquid_mixture_DIPPR9H,
        }
        for k, v in vals.items():
            print(f'Testing {k}')
            f = getattr(chemicals.thermal_conductivity, k)
            args = inspect.getfullargspec(f)[0]
            # Get the first element from each prop array
            kwargs = {i: np.atleast_1d(phase[argmap[i]])[0] for i in args}
            # If a dict was returned, take first element from each and put into list
            for item in list(kwargs.keys()):
                if hasattr(kwargs[item], 'keys'):
                    kwargs[item] = \
                        [kwargs[item][comp][0] for comp in phase.components.keys()]
            ref = f(**kwargs)
            val = v(phase)[0]
            assert_allclose(ref, val, rtol=1e-13)


if __name__ == '__main__':

    t = ThermalConductivityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()
