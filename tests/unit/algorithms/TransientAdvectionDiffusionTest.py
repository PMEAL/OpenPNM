import numpy as np
import scipy as sp
import openpnm as op


class TransientAdvectionDiffusionTest:

    def setup_class(self):
        np.random.seed(0)
        self.net = op.network.Cubic(shape=[4, 3, 1], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)

        self.phase = op.phases.GenericPhase(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['throat.diffusive_conductance'] = 1e-15
        self.phys['throat.hydraulic_conductance'] = 1e-15
        self.geo['pore.volume'] = 1e-27
        self.geo['throat.conduit_lengths.pore1'] = 0.1
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.1

    def test_transient_advection_diffusion(self):
        sf = op.algorithms.StokesFlow(network=self.net, phase=self.phase)
        sf.settings.update({'quantity': 'pore.pressure',
                            'conductance': 'throat.hydraulic_conductance'})
        sf.set_value_BC(pores=self.net.pores('right'), values=1)
        sf.set_value_BC(pores=self.net.pores('left'), values=0)
        sf.run()
        self.phase[sf.settings['quantity']] = sf[sf.settings['quantity']]

        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phys.add_model(propname='throat.ad_dif_conductance', model=mod,
                            s_scheme='powerlaw')
        self.phys.regenerate_models()

        ad = op.algorithms.TransientAdvectionDiffusion(network=self.net,
                                                       phase=self.phase)
        ad.settings.update({'phase': self.phase.name,
                            'quantity': 'pore.concentration',
                            'conductance':'throat.ad_dif_conductance',
                            'diffusive_conductance': 'throat.diffusive_conductance',
                            'hydraulic_conductance': 'throat.hydraulic_conductance',
                            'pressure': 'pore.pressure',
                            't_initial': 0,
                            't_final': 100,
                            't_step': 1,
                            't_output': 50,
                            't_tolerance': 1e-20,
                            't_precision': 12,
                            's_schem': 'implicit'})
        ad.set_IC(0)
        ad.set_value_BC(pores=self.net.pores('right'), values=2)
        ad.set_value_BC(pores=self.net.pores('left'), values=0)
        ad.run()

        x = [0., 0., 0.,
             0.89653, 0.89653, 0.89653,
             1.53924, 1.53924, 1.53924,
             2., 2., 2.]
        y = np.around(ad[ad.settings['quantity']], decimals=5)
        assert np.all(x == y)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = TransientAdvectionDiffusionTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
