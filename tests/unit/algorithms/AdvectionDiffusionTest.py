import openpnm as op
import scipy as sp
from numpy.testing import assert_allclose


class AdvectionDiffusionTest:

    def setup_class(self):
        sp.random.seed(0)
        self.net = op.network.Cubic(shape=[4, 3, 1], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['throat.conduit_lengths.pore1'] = 0.1
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.1

        self.phase = op.phases.GenericPhase(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['throat.diffusive_conductance'] = 1e-15
        self.phys['throat.hydraulic_conductance'] = 1e-15

        self.sf = op.algorithms.StokesFlow(network=self.net, phase=self.phase)
        self.sf.set_value_BC(pores=self.net.pores('back'), values=1)
        self.sf.set_value_BC(pores=self.net.pores('front'), values=0)
        self.sf.run()

        self.phase.update(self.sf.results())

        self.ad = op.algorithms.AdvectionDiffusion(network=self.net,
                                                   phase=self.phase)
        self.ad.settings.update({"cache_A": False, "cache_b": False})
        self.ad.set_value_BC(pores=self.net.pores('back'), values=2)
        self.ad.set_value_BC(pores=self.net.pores('front'), values=0)

    def test_conductance_automatically_gets_updated_after_StokesFlow(self):
        r"""
        This test ensures that ad_dif_conductance automatically gets updated
        when calling AdvectionDiffusion.run(), otherwise, a uniform pressure
        field will be used inside ad_dif_conductance method. Alternatively,
        one can manually call regenerate_models() after running StokesFlow,
        which was painful, hence ad_dif_conductance is now treated as an
        iterative_prop, so it automatically gets updated in the process.
        """
        # Add ad_dif pore-scale model prior to running StokesFlow
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phys.add_model(propname='throat.ad_dif_conductance',
                            model=mod, s_scheme='powerlaw')
        # Set up and run StokesFlow
        sf = op.algorithms.StokesFlow(network=self.net, phase=self.phase)
        sf.set_value_BC(pores=self.net.pores("back"), values=2)
        sf.set_value_BC(pores=self.net.pores("front"), values=0.0)
        sf.run()
        self.phase.update(sf.results())
        # Set up and run AdvectionDiffusion
        fd = op.algorithms.AdvectionDiffusion(network=self.phase, phase=self.phase)
        fd.set_value_BC(pores=self.net.pores("back"), values=1.0)
        fd.set_value_BC(pores=self.net.pores("front"), values=0.0)
        fd.run()
        # Verify that the new pressure field has been used
        c = fd["pore.concentration"]
        c_desired = [0, 0, 0, 0.56164952, 0.56164952, 0.56164952,
                     0.85096674, 0.85096674, 0.85096674, 1, 1, 1]
        assert_allclose(actual=c, desired=c_desired)
        self.setup_class()

    def test_powerlaw_advection_diffusion(self):
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phys.add_model(propname='throat.ad_dif_conductance_powerlaw',
                            model=mod, s_scheme='powerlaw')
        self.phys.regenerate_models()

        self.ad.setup(conductance='throat.ad_dif_conductance_powerlaw')
        self.ad.run()
        x = [0., 0., 0.,
             0.89653, 0.89653, 0.89653,
             1.53924, 1.53924, 1.53924,
             2., 2., 2.]
        y = sp.around(self.ad['pore.concentration'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_upwind_advection_diffusion(self):
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phys.add_model(propname='throat.ad_dif_conductance_upwind',
                            model=mod, s_scheme='upwind')
        self.phys.regenerate_models()
        self.ad.setup(conductance='throat.ad_dif_conductance_upwind')
        self.ad.run()
        x = [0., 0., 0.,
             0.86486, 0.86486, 0.86486,
             1.51351, 1.51351, 1.51351,
             2., 2., 2.]
        y = sp.around(self.ad['pore.concentration'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_hybrid_advection_diffusion(self):
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phys.add_model(propname='throat.ad_dif_conductance_hybrid',
                            model=mod, s_scheme='hybrid')
        self.phys.regenerate_models()

        self.ad.setup(conductance='throat.ad_dif_conductance_hybrid')
        self.ad.run()
        x = [0., 0., 0.,
             0.89908, 0.89908, 0.89908,
             1.54128, 1.54128, 1.54128,
             2., 2., 2.]
        y = sp.around(self.ad['pore.concentration'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_exponential_advection_diffusion(self):
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phys.add_model(propname='throat.ad_dif_conductance_exponential',
                            model=mod, s_scheme='exponential')
        self.phys.regenerate_models()

        self.ad.setup(conductance='throat.ad_dif_conductance_exponential')
        self.ad.run()
        x = [0., 0., 0.,
             0.89688173, 0.89688173, 0.89688173,
             1.53952557, 1.53952557, 1.53952557,
             2., 2., 2.]
        y = self.ad['pore.concentration']
        assert_allclose(actual=y, desired=x)

    def test_outflow_BC(self):
        for s_scheme in ['upwind', 'hybrid', 'powerlaw', 'exponential']:
            ad = op.algorithms.AdvectionDiffusion(network=self.net,
                                                  phase=self.phase)
            ad.setup(quantity='pore.concentration',
                     conductance='throat.ad_dif_conductance_'+s_scheme)

            ad.set_value_BC(pores=self.net.pores('back'), values=2)
            ad.set_outflow_BC(pores=self.net.pores('front'))
            ad.run()

            y = ad[ad.settings['quantity']].mean()
            assert_allclose(actual=y, desired=2.0)

    def test_rate(self):
        for s_scheme in ['upwind', 'hybrid', 'powerlaw', 'exponential']:
            ad = op.algorithms.AdvectionDiffusion(network=self.net,
                                                  phase=self.phase)
            ad.setup(quantity='pore.concentration',
                     conductance='throat.ad_dif_conductance_'+s_scheme,
                     s_scheme=s_scheme)

            ad.set_value_BC(pores=self.net.pores('back'), values=2)
            ad.set_value_BC(pores=self.net.pores('front'), values=0)
            ad.run()

            mdot_inlet = ad.rate(pores=self.net.pores("back"))[0]
            mdot_outlet = ad.rate(pores=self.net.pores("front"))[0]
            temp = sp.random.choice(self.net.pores(["back", "front"],
                                                   mode="not"),
                                    size=3, replace=False)
            mdot_internal = ad.rate(pores=temp)[0]
            # Ensure no mass is generated within the network
            assert_allclose(mdot_inlet - mdot_internal, mdot_inlet)
            # Ensure mass is not conserved
            assert_allclose(mdot_inlet, -mdot_outlet)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':
    t = AdvectionDiffusionTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
