import pytest
import numpy as np
import scipy as sp
import openpnm as op
from openpnm.phases import mixtures


class MixtureTest:
    def setup_class(self):
        ws = op.Workspace()
        ws.clear()
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.N2 = mixtures.species.gases.N2(network=self.net, name='pure_N2')
        self.O2 = mixtures.species.gases.O2(network=self.net, name='pure_O2')
        self.CO2 = mixtures.species.gases.CO2(network=self.net, name='pure_CO2')
        self.H2 = mixtures.species.gases.H2(network=self.net, name='pure_H2')
        self.air = mixtures.GenericMixture(network=self.net,
                                           components=[self.N2, self.O2,
                                                       self.H2, self.CO2],
                                           name='air_mixture')

    def test_set_mole_fraction(self):
        self.air.set_mole_fraction(self.N2, 0.790)
        self.air.set_mole_fraction(self.O2, 0.209)
        self.air.set_mole_fraction(self.CO2, 0.001)
        assert np.all(self.air['pore.mole_fraction.all'] == 1.0)

    def test_props(self):
        a = self.air.props(deep=False)
        b = self.air.props(deep=True)
        assert len(b) > len(a)

    def test_update_mole_fraction_with_molar_density(self):
        self.air.pop('pore.concentration.'+self.N2.name, None)
        self.air['pore.concentration.'+self.O2.name] = 0.5
        self.air['pore.concentration.'+self.CO2.name] = 0.0
        self.air['pore.concentration.'+self.H2.name] = 0.0
        self.air['pore.molar_density'] = 2.0
        self.air.update_mole_fractions(molar_density='pore.molar_density')
        assert np.all(self.air['pore.mole_fraction.all'] == 1.0)

    def test_update_mole_fraction_with_all_concentrations(self):
        self.air['pore.concentration.'+self.O2.name] = 1.5
        self.air['pore.concentration.'+self.N2.name] = 0.5
        self.air['pore.concentration.'+self.CO2.name] = 0.0
        self.air['pore.concentration.'+self.H2.name] = 0.0
        self.air.update_mole_fractions()
        assert np.all(self.air['pore.mole_fraction.all'] == 1.0)

    def test_interleave_data(self):
        r"""
        """
        self.air['pore.concentration.'+self.O2.name] = 1.0
        self.air['pore.concentration.'+self.N2.name] = 0.0
        self.air['pore.concentration.'+self.CO2.name] = 0.0
        self.air['pore.concentration.'+self.H2.name] = 0.0
        self.air.update_mole_fractions()
        MW = self.air['pore.molecular_weight'][0]
        assert MW == 0.032
        self.air['pore.concentration.'+self.O2.name] = 1.0
        self.air['pore.concentration.'+self.N2.name] = 0.5
        self.air['pore.concentration.'+self.CO2.name] = 0.3
        self.air['pore.concentration.'+self.H2.name] = 0.2
        self.air.update_mole_fractions()
        MW = self.air['pore.molecular_weight'][0]
        assert MW == 0.0298031

    def test_check_health(self):
        self.air.set_mole_fraction(self.N2, 0.790)
        self.air.set_mole_fraction(self.O2, 0.209)
        self.air.set_mole_fraction(self.CO2, 0.001)
        self.air.set_mole_fraction(self.H2, 0.000)
        h = self.air.check_mixture_health()
        assert h.health is True
        self.air.set_mole_fraction(self.CO2, 0.002)
        h = self.air.check_mixture_health()
        assert h.health is False
        self.air.set_mole_fraction(self.CO2, 0.000)
        h = self.air.check_mixture_health()
        assert h.health is False
        assert len(h['mole_fraction_too_low']) == self.air.Np
        self.air.set_mole_fraction(self.CO2, 0.001)
        self.air['pore.mole_fraction.'+self.CO2.name][0] = 0.0
        h = self.air.check_mixture_health()
        assert h.health is False
        assert len(h['mole_fraction_too_low']) == 1

    def test_getitem(self):
        d = self.air['pore.mole_fraction']
        set_a = set(['pore.mole_fraction.pure_N2',
                     'pore.mole_fraction.pure_O2',
                     'pore.mole_fraction.pure_H2',
                     'pore.mole_fraction.pure_CO2',
                     'pore.mole_fraction.all'])
        assert set_a.difference(set(d.keys())) == set()


if __name__ == '__main__':

    t = MixtureTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
