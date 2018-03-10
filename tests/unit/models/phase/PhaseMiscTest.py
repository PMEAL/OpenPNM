import OpenPNM
import scipy as sp
import OpenPNM.Phases.models as phm


class PhaseMiscTest:

    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.air = OpenPNM.Phases.Air(network=self.net)

    def test_linear(self):
        catch = self.air.pop('pore.test', None)
        self.air['pore.constant'] = 2
        mod = phm.misc.linear
        self.air.models.add(model=mod,
                            propname='pore.test',
                            poreprop='pore.constant',
                            m=2,
                            b=1)

        assert sp.all(self.air['pore.test'] == 5)

    def test_polynomial_2terms(self):
        catch = self.air.pop('pore.test', None)
        self.air['pore.constant'] = 2
        mod = phm.misc.polynomial
        self.air.models.add(model=mod,
                            propname='pore.test',
                            poreprop='pore.constant',
                            a=[1, 2])

        assert sp.all(self.air['pore.test'] == 5)

    def test_polynomial_3terms(self):
        catch = self.air.pop('pore.test', None)
        self.air['pore.constant'] = 2
        mod = phm.misc.polynomial
        self.air.models.add(model=mod,
                            propname='pore.test',
                            poreprop='pore.constant',
                            a=[1, 2, 3])

        assert sp.all(self.air['pore.test'] == 17)

    def test_ideal_mixture(self):
        pass

    def test_mixture_value(self):
        pass


if __name__ == '__main__':
    a = PhaseMiscTest()
    a.setup_class()
    self = a
