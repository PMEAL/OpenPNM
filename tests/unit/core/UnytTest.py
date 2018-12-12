import openpnm as op
import scipy as sp
import pytest
import unyt


class UnytTest:

    def setup_class(self):
        self.ws = op.Workspace()
        self.ws.clear()
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.geo1 = op.geometry.StickAndBall(network=self.net,
                                             pores=sp.arange(0, 75),
                                             throats=sp.arange(0, 150))
        self.geo2 = op.geometry.StickAndBall(network=self.net,
                                             pores=sp.arange(75, self.net.Np),
                                             throats=sp.arange(150, self.net.Nt))
        self.phase = op.phases.Air(network=self.net)
        self.phys1 = op.physics.Standard(network=self.net, phase=self.phase,
                                         geometry=self.geo1)
        self.phys2 = op.physics.Standard(network=self.net, phase=self.phase,
                                         geometry=self.geo2)

    def test_set_and_get_on_network(self):
        self.net['pore.test_prop'] = sp.ones(self.net.Np) * unyt.m
        t = self.net['pore.test_prop']
        assert t.units == unyt.m

    def test_set_and_get_on_phase(self):
        self.phase['pore.test_prop2'] = sp.ones(self.phase.Np) * unyt.m
        t = self.phase['pore.test_prop2']
        assert t.units == unyt.m

    def test_set_on_net_get_on_geo1_and_geo2(self):
        self.net['pore.test_prop'] = sp.ones(self.net.Np) * unyt.m
        t = self.geo1['pore.test_prop']
        assert t.units == unyt.m
        t = self.geo2['pore.test_prop']
        assert t.units == unyt.m

    def test_set_on_phase_get_on_phys1_and_phys2(self):
        self.phase['pore.test_prop2'] = sp.ones(self.phase.Np) * unyt.m
        t = self.phys1['pore.test_prop2']
        assert t.units == unyt.m
        t = self.phys2['pore.test_prop2']
        assert t.units == unyt.m

    def test_set_on_geo1_and_geo2_get_on_net(self):
        self.geo1['pore.test_prop3'] = sp.ones(self.geo1.Np) * unyt.m
        self.geo2['pore.test_prop3'] = sp.ones(self.geo2.Np) * unyt.m
        t = self.net['pore.test_prop3']
        assert t.units == unyt.m

    def test_set_on_geo1_not_geo2_get_on_net(self):
        self.geo1['pore.test_prop4'] = sp.ones(self.geo1.Np) * unyt.m
        t = self.net['pore.test_prop4']
        assert t.units == unyt.m

    def test_set_on_phys1_and_phys2_get_on_phase(self):
        self.phys1['pore.test_prop5'] = sp.ones(self.phys1.Np) * unyt.m
        self.phys2['pore.test_prop5'] = sp.ones(self.phys2.Np) * unyt.m
        t = self.phase['pore.test_prop5']
        assert t.units == unyt.m

    def test_set_on_phys1_not_phys2_get_on_phase(self):
        self.phys1['pore.test_prop6'] = sp.ones(self.phys1.Np) * unyt.m
        t = self.phase['pore.test_prop6']
        assert t.units == unyt.m

    def test_interleave_compatible_units(self):
        self.phys1['pore.test_prop5'] = sp.ones(self.phys1.Np) * unyt.m
        self.phys2['pore.test_prop5'] = sp.ones(self.phys2.Np) * unyt.cm
        t = self.phase['pore.test_prop5']
        assert t.units == unyt.m

    def test_interleave_incompatible_units(self):
        self.phys1['pore.test_prop5'] = sp.ones(self.phys1.Np) * unyt.m
        self.phys2['pore.test_prop5'] = sp.ones(self.phys2.Np) * unyt.kg
        with pytest.raises(Exception):
            self.phase['pore.test_prop5']

    def test_add_units_to_models(self):
        self.geo1['pore.diameter'] *= unyt.m
        self.geo1.regenerate_models(propnames=['pore.volume'])
        t = self.geo1['pore.volume']
        assert t.units == unyt.m**3


if __name__ == '__main__':

    t = UnytTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
