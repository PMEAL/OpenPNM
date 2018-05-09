import openpnm as op
import scipy as sp
import pytest


class TransientCranknicolsonReactiveTransportTest:

    def setup_class(self):
        sp.random.seed(0)
        self.net = op.network.Cubic(shape=[3, 3, 1], spacing=1e-4)
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)

        self.phase = op.phases.Water(network=self.net)
        self.phase['throat.viscosity'] = self.phase['pore.viscosity'][0]

        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['pore.A'] = 1e-10
        self.phys['pore.k'] = 2
        mod1 = op.models.physics.hydraulic_conductance.hagen_poiseuille
        mod2 = op.models.physics.generic_source_term.standard_kinetics
        self.phys.add_model(propname='throat.conductance',
                            model=mod1,
                            viscosity='throat.viscosity',
                            regen_mode='normal')
        self.phys.add_model(propname='pore.reaction',
                            model=mod2,
                            prefactor='pore.A',
                            exponent='pore.k',
                            quantity='pore.pressure',
                            regen_mode='normal')
        self.s = {'conductance': 'throat.conductance',
                  'quantity': 'pore.pressure'}

    def test_transient_cranknicolson_reactive_transport(self):
        alg = op.algorithms.TransientReactiveTransport(network=self.net,
                                                       phase=self.phase,
                                                       settings=self.s)
        alg.settings.update({'t_scheme': 'cranknicolson', 't_step': 0.1,
                             't_tolerance': 1e-07, 'r_tolerance': 1e-06})
        alg.set_IC(0)
        alg.set_value_BC(pores=self.net.pores('left'), values=2)
        alg.set_source(propname='pore.reaction', pores=self.net.pores('right'))
        alg.run()
        x = [2, 8.0339e-01, 4.4e-04,
             2, 7.3107e-01, 2.1e-04,
             2, 1.3544e-01, 4e-04]
        y = sp.around(alg[alg.settings['quantity']], decimals=5)
        assert sp.all(x == y)

    def teardown_class(self):
        ws = op.core.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = TransientCranknicolsonReactiveTransportTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
