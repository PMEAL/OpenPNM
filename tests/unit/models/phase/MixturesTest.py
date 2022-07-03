import openpnm as op
import numpy as np
import openpnm.models as mods


class MixturesTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
    #     self.N2 = mixtures.species.gases.N2(network=self.net)
    #     self.O2 = mixtures.species.gases.O2(network=self.net)
    #     self.H2O = mixtures.species.liquids.H2O(network=self.net)
    #     self.ha = mixtures.GenericMixture(network=self.net,
    #                                       components=[self.N2,
    #                                                   self.O2,
    #                                                   self.H2O])

    # def test_wilke_fuller_diffusivity(self):
    #     self.ha.set_concentration(component=self.N2, values=0.7)
    #     self.ha.set_concentration(component=self.O2, values=0.2)
    #     self.ha.set_concentration(component=self.H2O, values=0.1)
    #     self.ha.update_mole_fractions()
    #     self.ha.add_model(propname='pore.diffusivity',
    #                       model=mods.phase.mixtures.wilke_fuller_diffusivity)
    #     assert 'pore.diffusivity.' + self.O2.name in self.ha.keys()
    #     assert 'pore.diffusivity.' + self.N2.name in self.ha.keys()
    #     assert 'pore.diffusivity.' + self.H2O.name in self.ha.keys()
    #     Da = self.ha['pore.diffusivity.' + self.O2.name][0]
    #     Db = self.ha['pore.diffusivity.' + self.N2.name][0]
    #     Dc = self.ha['pore.diffusivity.' + self.H2O.name][0]
    #     assert Da != Db
    #     assert Db != Dc

    #     self.ha.set_concentration(component=self.N2, values=0.8)
    #     self.ha.set_concentration(component=self.O2, values=0.2)
    #     self.ha.set_concentration(component=self.H2O, values=0.0)
    #     self.ha.update_mole_fractions()
    #     self.ha.regenerate_models()
    #     Da = self.ha['pore.diffusivity.' + self.O2.name][0]
    #     Db = self.ha['pore.diffusivity.' + self.N2.name][0]
    #     Dc = self.ha['pore.diffusivity.' + self.H2O.name][0]
    #     assert np.allclose(Da, Db)
    #     assert Db != Dc

    # def test_mole_weighted_average(self):
    #     self.ha.set_concentration(component=self.N2, values=0.5)
    #     self.ha.set_concentration(component=self.O2, values=0.5)
    #     self.ha.set_concentration(component=self.H2O, values=0.0)
    #     self.N2['pore.test'] = 3
    #     self.O2['pore.test'] = 3
    #     self.ha.add_model(propname='pore.blah',
    #                       model=mods.phase.mixtures.mole_weighted_average,
    #                       prop='pore.test')
    #     assert 'pore.blah' not in self.ha.keys()
    #     self.H2O['pore.test'] = 3
    #     self.ha.add_model(propname='pore.blah',
    #                       model=mods.phase.mixtures.mole_weighted_average,
    #                       prop='pore.test')
    #     assert 'pore.blah' in self.ha.keys()


if __name__ == '__main__':

    t = MixturesTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
