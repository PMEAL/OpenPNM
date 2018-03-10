import openpnm as op
import scipy as sp
import openpnm.models.phase as phm


class PhaseMiscTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.air = op.phases.Air(network=self.net)

    def test_linear(self):
        catch = self.air.pop('pore.test', None)
        self.air['pore.constant'] = 2
        self.air.add_model(model=phm.misc.linear,
                           propname='pore.test',
                           prop='pore.constant',
                           m=2, b=1)
        self.air.regenerate_models()
        assert sp.all(self.air['pore.test'] == 5)

    def test_polynomial_2terms(self):
        catch = self.air.pop('pore.test', None)
        self.air['pore.constant'] = 2
        self.air.add_model(model=phm.misc.polynomial,
                           propname='pore.test',
                           prop='pore.constant',
                           a=[1, 2])
        self.air.regenerate_models()
        assert sp.all(self.air['pore.test'] == 5)

    def test_polynomial_3terms(self):
        catch = self.air.pop('pore.test', None)
        self.air['pore.constant'] = 2
        self.air.add_model(model=phm.misc.polynomial,
                           propname='pore.test',
                           prop='pore.constant',
                           a=[1, 2, 3])
        self.air.regenerate_models()
        assert sp.all(self.air['pore.test'] == 17)

    def test_ideal_mixture(self):
        pass

    def test_mixture_value(self):
        pass


if __name__ == '__main__':

    t = PhaseMiscTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
