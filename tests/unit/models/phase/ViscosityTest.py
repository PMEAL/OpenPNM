import openpnm as op
import chemicals
from thermo import Chemical
from numpy.testing import assert_array_almost_equal, assert_allclose
from openpnm.utils import get_mixture_model_args


class ViscosityTest:

    def test_water_correlation(self):
        pn = op.network.Demo()
        water = op.phase.Water(network=pn)
        temp = Chemical('h2o')
        assert_allclose(temp.mu, water['pore.viscosity'][0], rtol=0.03)

    def test_pure_gas(self):
        pn = op.network.Demo()
        gas = op.phase.Species(network=pn, species='ch4')
        gas.add_model(propname='pore.viscosity',
                      model=op.models.phase.viscosity.gas_pure)
        temp = chemicals.viscosity.viscosity_gas_Gharagheizi(
            T=gas['pore.temperature'][0],
            Tc=gas['param.critical_temperature'],
            Pc=gas['param.critical_pressure'],
            MW=gas['param.molecular_weight'],
        )
        assert_allclose(temp, gas['pore.viscosity'][0], rtol=1e-10)

    def test_pure_liquid(self):
        pn = op.network.Demo()
        liq = op.phase.Species(network=pn, species='h2o')
        liq.add_model(propname='pore.viscosity',
                      model=op.models.phase.viscosity.liquid_pure)
        temp = chemicals.viscosity.Letsou_Stiel(
            T=liq['pore.temperature'][0],
            Tc=liq['param.critical_temperature'],
            Pc=liq['param.critical_pressure'],
            MW=liq['param.molecular_weight'],
            omega=liq['param.acentric_factor'],
        )
        assert_allclose(temp, liq['pore.viscosity'][0], rtol=1e-10)

    def test_gas_mixture(self):
        pn = op.network.Demo()
        gas1 = op.phase.Species(network=pn, species='ch4')
        gas2 = op.phase.Species(network=pn, species='co2')
        gas1.add_model(propname='pore.viscosity',
                       model=op.models.phase.viscosity.gas_pure)
        gas2.add_model(propname='pore.viscosity',
                       model=op.models.phase.viscosity.gas_pure)
        mix = op.phase.GasMixture(network=pn, components=[gas1, gas2])
        mix.y(gas1, 0.5)
        mix.y(gas2, 0.5)
        mix.add_model(propname='pore.viscosity',
                      model=op.models.phase.viscosity.gas_mixture)
        args = get_mixture_model_args(mix, composition='zs',
                                      args={'mus': 'pore.viscosity',
                                            'MWs': 'param.molecular_weight'})
        temp = chemicals.viscosity.Herning_Zipperer(**args)
        assert_allclose(temp, mix['pore.viscosity'][0], rtol=1e-10)

    def test_liquid_mixture(self):
        pn = op.network.Demo()
        liq1 = op.phase.Species(network=pn, species='h2o')
        liq2 = op.phase.Species(network=pn, species='etoh')
        liq1.add_model(propname='pore.viscosity',
                       model=op.models.phase.viscosity.liquid_pure)
        liq2.add_model(propname='pore.viscosity',
                       model=op.models.phase.viscosity.liquid_pure)
        mix = op.phase.LiquidMixture(network=pn, components=[liq1, liq2])
        mix.x(liq1, 0.5)
        mix.x(liq2, 0.5)
        mix.add_model(propname='pore.viscosity',
                      model=op.models.phase.viscosity.liquid_mixture)
        args = get_mixture_model_args(mix, composition='fracs',
                                      args={'props': 'pore.viscosity'})
        temp = chemicals.utils.mixing_simple(**args)
        assert_allclose(temp, mix['pore.viscosity'][0], rtol=0.02)

    def test_generic_chemicals_for_pure_gas_viscosity(self):
        mods = [
            'viscosity_gas_Gharagheizi',
            'Yoon_Thodos',
            'Stiel_Thodos',
            # 'Lucas_gas',  # This one does not work
        ]
        pn = op.network.Demo()
        o2 = op.phase.Species(network=pn, species='oxygen')
        mu = []
        for m in mods:
            f = getattr(chemicals.viscosity, m)
            mu.append(op.models.phase.chemicals_wrapper(target=o2, f=f).mean())
        assert_array_almost_equal(mu, 2.05e-5, decimal=6)

    def test_generic_chemicals_for_pure_liq_viscosity(self):
        mods = [
            # The following 3 require constants
            # 'Viswanath_Natarajan_3',
            # 'Viswanath_Natarajan_2',
            # 'Viswanath_Natarajan_2_exponential',
            'Letsou_Stiel',
            # 'Przedziecki_Sridhar',  # Requires molar volume at temperature
            # 'Lucas',  # Requires saturation pressure at temperature
        ]
        pn = op.network.Demo()
        h2o = op.phase.Species(network=pn, species='water')
        mu = []
        for m in mods:
            f = getattr(chemicals.viscosity, m)
            mu.append(op.models.phase.chemicals_wrapper(target=h2o, f=f).mean())
        assert_allclose(mu, 0.000621, rtol=0.02)


if __name__ == '__main__':

    t = ViscosityTest()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
