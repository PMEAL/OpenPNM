import numpy as np
import pytest

import openpnm as op


class MixtureTest:
    def setup_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_init_binary_gas(self):
        net = op.network.Demo()
        N2 = op.phase.Species(network=net, species='n2', name='pure_N2')
        O2 = op.phase.Species(network=net, species='o2', name='pure_O2')
        air = op.phase.BinaryGas(network=net, components=[N2, O2], name='air')

    def test_liquid_mixture_two_comps(self):
        net = op.network.Demo()
        A = op.phase.StandardLiquid(network=net, species='ethanol')
        B = op.phase.StandardLiquid(network=net, species='h2o')
        vodka = op.phase.LiquidMixture(network=net, components=[A, B])
        vodka.x(A.name, 0.4)
        vodka.x(B.name, 0.6)
        assert len(vodka.components) == 2

    def test_standard_liquid_mixture(self):
        net = op.network.Demo()
        A = op.phase.StandardLiquid(network=net, species='ethanol')
        B = op.phase.StandardLiquid(network=net, species='h2o')
        vodka = op.phase.StandardLiquidMixture(network=net, components=[A, B])
        vodka.x(A.name, 0.4)
        vodka.x(B.name, 0.6)
        # Ensure models are NOT run during init (no point without mol fracs)
        m = vodka.models.keys()
        for item in m:
            assert item.split('@')[0] not in vodka.keys()
        # Make sure models run without complaining
        vodka.regenerate_models()
        # Make sure all models were actually run
        m = vodka.models.keys()
        for item in m:
            assert item.split('@')[0] in vodka.keys()

        vodka.add_model(propname='pore.diffusivity',
                        model=op.models.phase.diffusivity.liquid_mixture_tc)

    def test_standard_gas_mixture(self):
        net = op.network.Demo()
        A = op.phase.StandardGas(network=net, species='o2')
        B = op.phase.StandardGas(network=net, species='n2')
        air = op.phase.StandardGasMixture(network=net, components=[A, B])
        air.y(A.name, 0.21)
        air.y(B.name, 0.79)
        # Ensure models are NOT run during init (no point without mol fracs)
        m = air.models.keys()
        for item in m:
            assert item.split('@')[0] not in air.keys()
        # Make sure models run without complaining
        air.regenerate_models()
        # Make sure all models were actually run
        m = air.models.keys()
        for item in m:
            assert item.split('@')[0] in air.keys()

        air.add_model(propname='pore.diffusivity',
                      model=op.models.phase.diffusivity.gas_mixture_ce)

    def test_add_and_remove_component_method(self):
        net = op.network.Demo()
        o2 = op.phase.Species(network=net, species='o2', name='pure_O2')
        n2 = op.phase.Species(network=net, species='n2', name='pure_N2')
        air = op.phase.GasMixture(network=net, components=[n2, o2])
        air.remove_comp(n2)
        assert len(air.components) == 1
        air.remove_comp(o2)
        assert len(air.components) == 0

    def test_check_health(self):
        net = op.network.Demo()
        o2 = op.phase.Species(network=net, species='o2', name='pure_O2')
        n2 = op.phase.Species(network=net, species='n2', name='pure_N2')
        air = op.phase.GasMixture(network=net, components=[n2, o2])
        air['pore.mole_fraction.pure_N2'] = 0.79
        air['pore.mole_fraction.pure_O2'] = 0.21
        h = air.check_mixture_health()
        assert h.health is True
        air['pore.mole_fraction.pure_O2'] = 0.25
        h = air.check_mixture_health()
        assert h.health is False

    def test_getitem(self):
        net = op.network.Demo()
        o2 = op.phase.Species(network=net, species='o2', name='pure_O2')
        n2 = op.phase.Species(network=net, species='n2', name='pure_N2')
        air = op.phase.GasMixture(network=net, components=[n2, o2])
        d = air['pore.mole_fraction']
        set_a = set(['pure_N2', 'pure_O2'])
        assert set_a.difference(set(d.keys())) == set()

    def test_getitem_with_species(self):
        net = op.network.Demo()
        o2 = op.phase.StandardGas(network=net, species='o2', name='pure_O2')
        n2 = op.phase.StandardGas(network=net, species='n2', name='pure_N2')
        air = op.phase.GasMixture(network=net, components=[n2, o2])
        mu1 = air['pore.viscosity.pure_N2']
        mu2 = air['pore.viscosity.pure_O2']
        assert not np.all(mu1 == mu2)

    def test_get_comp(self):
        net = op.network.Demo()
        o2 = op.phase.StandardGas(network=net, species='o2', name='pure_O2')
        n2 = op.phase.StandardGas(network=net, species='n2', name='pure_N2')
        air = op.phase.GasMixture(network=net, components=[n2, o2])
        del o2['pore.viscosity']
        _ = air['pore.viscosity.pure_N2']
        with pytest.raises(KeyError):
            _ = air['pore.viscosity.pure_O2']

    def test_get_mix_vals(self):
        net = op.network.Demo()
        o2 = op.phase.StandardGas(network=net, species='o2', name='pure_O2')
        n2 = op.phase.StandardGas(network=net, species='n2', name='pure_N2')
        air = op.phase.GasMixture(network=net, components=[n2, o2])
        mu1 = air.get_mix_vals('pore.viscosity')
        # Mole fraction weighting won't work without mole fractions
        assert np.all(np.isnan(mu1))
        air.y(o2.name, 0.5)
        air.y(n2.name, 0.5)
        mu2 = air.get_mix_vals('pore.viscosity')
        assert np.all(np.isreal(mu2))
        # Assert each method works and returns something different
        mu3 = air.get_mix_vals('pore.viscosity', mode='logarithmic')
        mu4 = air.get_mix_vals('pore.viscosity', mode='power', power=2)
        assert not np.all(mu2 == mu3)
        assert not np.all(mu3 == mu4)
        assert not np.all(mu2 == mu4)

    def test_regenerate_components(self):
        net = op.network.Demo()
        o2 = op.phase.StandardGas(network=net, species='o2', name='pure_O2')
        n2 = op.phase.StandardGas(network=net, species='n2', name='pure_N2')
        air = op.phase.StandardGasMixture(network=net, components=[n2, o2])
        air.y(o2.name, y=0.2)
        air.y(n2.name, y=0.8)
        air.regenerate_models()
        mu1 = air['pore.viscosity'].mean()
        o2['pore.temperature'] = 333.3
        air.regenerate_models()
        mu2 = air['pore.viscosity'].mean()
        assert mu2 > mu1
        o2['pore.temperature'] = 298.0
        o2.regenerate_models()
        mu3 = air['pore.viscosity'].mean()
        assert mu2 == mu3
        air.regenerate_models()
        mu4 = air['pore.viscosity'].mean()
        assert mu4 == mu1


if __name__ == '__main__':

    t = MixtureTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()
