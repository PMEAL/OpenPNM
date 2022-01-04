import openpnm as op


class MixtureTest:
    def setup_class(self):
        ws = op.Workspace()
        ws.clear()
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.N2 = op.phase.GasByName(network=self.net, species='n2', name='pure_N2')
        self.O2 = op.phase.GasByName(network=self.net, species='o2', name='pure_O2')
        self.air = op.phase.GasMixture(network=self.net,
                                       components=[self.N2, self.O2],
                                       name='air_mixture')

    def test_props(self):
        a = self.air.props(deep=False)
        b = self.air.props(deep=True)
        assert len(b) > len(a)


    def test_check_health(self):
        self.air['pore.mole_fraction.pure_N2'] = 0.79
        self.air['pore.mole_fraction.pure_O2'] = 0.21
        h = self.air.check_mixture_health()
        assert h.health is True
        self.air['pore.mole_fraction.pure_O2'] = 0.25
        h = self.air.check_mixture_health()
        assert h.health is False


    def test_getitem(self):
        d = self.air['pore.mole_fraction']
        set_a = set(['pore.mole_fraction.pure_N2',
                     'pore.mole_fraction.pure_O2',
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
