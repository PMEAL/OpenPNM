import openpnm as op
import scipy as sp
import pytest


class TransientImplicitReactiveTransportTest:

    def setup_class(self):
        sp.random.seed(0)
        self.net = op.network.Cubic(shape=[3, 3, 1], spacing=1e-6)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.volume'] = 1e-15
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['pore.A'] = 1e-7
        self.phys['pore.k'] = 2
        self.phys['throat.diffusive_conductance'] = 1e-12
        mod = op.models.physics.generic_source_term.standard_kinetics
        self.phys.add_model(propname='pore.reaction',
                            model=mod,
                            prefactor='pore.A',
                            exponent='pore.k',
                            quantity='pore.concentration',
                            regen_mode='normal')
        self.s = {'conductance': 'throat.diffusive_conductance',
                  'quantity': 'pore.concentration'}

    def test_transient_implicit_reactive_transport(self):
        alg = op.algorithms.TransientReactiveTransport(network=self.net,
                                                       phase=self.phase,
                                                       settings=self.s)
        alg.settings.update({'t_scheme': 'implicit', 't_step': 0.1,
                             't_tolerance': 1e-07, 'r_tolerance': 1e-06})
        alg.set_IC(0)
        alg.set_value_BC(pores=self.net.pores('left'), values=2)
        alg.set_source(propname='pore.reaction', pores=self.net.pores('right'))
        alg.run()
        x = [2, 1.00158, 0.00316,
             2, 1.00158, 0.00316,
             2, 1.00158, 0.00316]
        y = sp.around(alg[alg.settings['quantity']], decimals=5)
        assert sp.all(x == y)

    def teardown_class(self):
        ws = op.core.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = TransientImplicitReactiveTransportTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
