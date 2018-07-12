import openpnm as op
import scipy as sp
import pytest


class GenericTransportTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[9, 9, 9])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['throat.diffusive_conductance'] = 1e-15
        self.phys['pore.A'] = -1e-10
        self.phys['pore.k'] = 2
        mod = op.models.physics.generic_source_term.standard_kinetics
        self.phys.add_model(propname='pore.reaction',
                            model=mod,
                            prefactor='pore.A',
                            exponent='pore.k',
                            quantity='pore.concentration',
                            regen_mode='normal')

    def test_one_value_one_source(self):
        rt = op.algorithms.ReactiveTransport(network=self.net,
                                             phase=self.phase)
        rt.settings.update({'conductance': 'throat.diffusive_conductance',
                            'quantity': 'pore.concentration'})
        rt.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        rt.set_value_BC(pores=self.net.pores('top'), values=1.0)
        rt.run()
        x = [0.0011, 0.1260, 0.2508, 0.3757, 0.5006, 0.6254, 0.7503, 0.8751, 1.0]
        y = sp.unique(sp.around(rt['pore.concentration'], decimals=4))
        assert sp.all(x == y)

    def test_source_not_override_BCs(self):
        rt = op.algorithms.ReactiveTransport(network=self.net,
                                             phase=self.phase)
        rt.settings.update({'conductance': 'throat.diffusive_conductance',
                            'quantity': 'pore.concentration'})
        rt.set_value_BC(pores=self.net.pores('left'), values=1.0)
        rt.set_value_BC(pores=self.net.pores('right'), values=0.5)
        rt.set_source(pores=self.net.Ps, propname='pore.reaction')
        rt.run()
        yleft = rt['pore.concentration'][self.net.pores('left')].mean()
        yright = rt['pore.concentration'][self.net.pores('right')].mean()
        assert yleft == 1.0
        assert yright == 0.5

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = GenericTransportTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
