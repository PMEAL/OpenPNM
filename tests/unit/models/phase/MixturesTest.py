import openpnm as op
from openpnm.phases import mixtures
import openpnm.models as mods


class MixturesTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])

    def test_wilke_fuller_diffusivity(self):
        N2 = mixtures.species.gases.N2(network=self.net)
        O2 = mixtures.species.gases.O2(network=self.net)
        H2O = mixtures.species.liquids.H2O(network=self.net)
        ha = mixtures.GenericMixture(network=self.net,
                                     components=[N2, O2, H2O])
        ha.set_concentration(component=N2, values=0.7)
        ha.set_concentration(component=O2, values=0.2)
        ha.set_concentration(component=H2O, values=0.1)
        ha.update_mole_fractions()
        ha.add_model(propname='pore.diffusivity',
                     model=mods.phases.mixtures.wilke_fuller_diffusivity)
        assert 'pore.diffusivity.' + O2.name in ha.keys()
        assert 'pore.diffusivity.' + N2.name in ha.keys()
        assert 'pore.diffusivity.' + H2O.name in ha.keys()


if __name__ == '__main__':

    t = MixturesTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
