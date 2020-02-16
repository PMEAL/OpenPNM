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
        self.phys['pore.A'] = -1e-7
        self.phys['pore.k'] = 2
        self.phys['throat.diffusive_conductance'] = 1e-12
        mod = op.models.physics.generic_source_term.standard_kinetics
        self.phys.add_model(propname='pore.reaction',
                            model=mod,
                            prefactor='pore.A',
                            exponent='pore.k',
                            X='pore.concentration',
                            regen_mode='deferred')
        self.settings = {'conductance': 'throat.diffusive_conductance',
                         'quantity': 'pore.concentration'}

    def test_transient_implicit_reactive_transport(self):
        alg = op.algorithms.TransientReactiveTransport(network=self.net,
                                                       phase=self.phase,
                                                       settings=self.settings)
        alg.setup(quantity='pore.concentration',
                  conductance='throat.diffusive_conductance',
                  t_initial=1, t_final=1000, t_step=0.1, t_output=100,
                  t_tolerance=1e-07, t_precision=14, t_scheme='implicit')
        alg.settings.update({'rxn_tolerance': 1e-06})
        alg.set_IC(0)
        alg.set_value_BC(pores=self.net.pores('left'), values=2)
        alg.set_source(propname='pore.reaction', pores=self.net.pores('right'))
        alg.run()
        x = [2., 1.00158, 0.00316,
             2., 1.00158, 0.00316,
             2., 1.00158, 0.00316]
        y = sp.around(alg[alg.settings['quantity']], decimals=5)
        assert sp.all(x == y)

    def test_transient_cranknicolson_reactive_transport(self):
        alg = op.algorithms.TransientReactiveTransport(network=self.net,
                                                       phase=self.phase,
                                                       settings=self.settings)
        alg.setup(t_initial=1, t_final=11, t_precision=14)
        alg.settings.update({'t_scheme': 'cranknicolson', 't_step': 0.1,
                             't_tolerance': 1e-07, 'rxn_tolerance': 1e-06})
        alg.set_value_BC(pores=self.net.pores('left'), values=2)
        alg.set_source(propname='pore.reaction', pores=self.net.pores('right'))
        alg.run()
        x = ([2., 1.02351, 0.04428,
              2., 1.02351, 0.04428,
              2., 1.02351, 0.04428])
        y = sp.around(alg[alg.settings['quantity']], decimals=5)
        assert sp.all(x == y)

    def test_transient_reactive_transport_output(self):
        alg = op.algorithms.TransientReactiveTransport(network=self.net,
                                                       phase=self.phase,
                                                       settings=self.settings)
        alg.setup(t_initial=2, t_final=12, t_precision=10)
        alg.settings.update({'t_scheme': 'cranknicolson', 't_step': 0.1,
                             't_tolerance': 1e-07, 'rxn_tolerance': 1e-06,
                             't_output': sp.arange(4, 7, 1)})
        alg.set_value_BC(pores=self.net.pores('left'), values=2)
        alg.set_source(propname='pore.reaction', pores=self.net.pores('right'))
        alg.run()
        times = ['pore.concentration@2', 'pore.concentration@4',
                 'pore.concentration@5', 'pore.concentration@6',
                 'pore.concentration@12']
        assert (set(times).issubset(set(alg.keys())))

    def test_transient_reactive_transport_results(self):
        alg = op.algorithms.TransientReactiveTransport(network=self.net,
                                                       phase=self.phase,
                                                       settings=self.settings)
        alg.setup(t_initial=2, t_final=12, t_precision=10)
        alg.settings.update({'t_scheme': 'cranknicolson', 't_step': 0.1,
                             't_tolerance': 1e-07, 'rxn_tolerance': 1e-06,
                             't_output': sp.arange(2, 13, 1)})
        alg.set_value_BC(pores=self.net.pores('left'), values=2)
        alg.set_source(propname='pore.reaction', pores=self.net.pores('right'))
        alg.run()
        times_1 = ['pore.concentration@5', 'pore.concentration@11']
        results_times_1 = alg.results(times=[5, 11]).keys()
        times_2 = ['pore.concentration@7']
        results_times_2 = alg.results(times=7).keys()
        times_3 = ['pore.concentration@2', 'pore.concentration@3',
                   'pore.concentration@4', 'pore.concentration@5',
                   'pore.concentration@6', 'pore.concentration@7',
                   'pore.concentration@8', 'pore.concentration@9',
                   'pore.concentration@10', 'pore.concentration@11',
                   'pore.concentration@12']
        results_times_3 = alg.results(steps=None).keys()
        assert (set(times_1).issubset(set(results_times_1))
                and set(times_2).issubset(set(results_times_2))
                and set(times_3).issubset(set(results_times_3)))

    def test_transient_steady_mode_reactive_transport(self):
        alg = op.algorithms.TransientReactiveTransport(network=self.net,
                                                       phase=self.phase,
                                                       settings=self.settings)
        alg.settings.update({'t_scheme': 'steady', 'rxn_tolerance': 1e-06})
        alg.set_IC(0)
        alg.set_value_BC(pores=self.net.pores('left'), values=2)
        alg.set_source(propname='pore.reaction', pores=self.net.pores('right'))
        alg.run()
        x = [2, 1.00158, 0.00316,
             2, 1.00158, 0.00316,
             2, 1.00158, 0.00316]
        y = sp.around(alg[alg.settings['quantity']], decimals=5)
        assert sp.all(x == y)

    def test_adding_bc_over_sources(self):
        alg = op.algorithms.TransientReactiveTransport(network=self.net,
                                                       phase=self.phase,
                                                       settings=self.settings)
        alg.set_source(propname='pore.reaction', pores=self.net.pores('right'))
        with pytest.raises(Exception):
            alg.set_value_BC(pores=self.net.pores('right'), values=2)

    def test_adding_sources_over_bc(self):
        alg = op.algorithms.TransientReactiveTransport(network=self.net,
                                                       phase=self.phase,
                                                       settings=self.settings)
        alg.set_value_BC(pores=self.net.pores('right'), values=2)
        with pytest.raises(Exception):
            alg.set_source(propname='pore.reaction',
                           pores=self.net.pores('right'))

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = TransientImplicitReactiveTransportTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
