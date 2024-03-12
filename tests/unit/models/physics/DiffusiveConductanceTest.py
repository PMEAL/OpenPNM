import numpy as np
import openpnm as op
from numpy.testing import assert_allclose


class DiffusiveConductanceTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.net['pore.diameter'] = 1.0
        self.net['throat.diameter'] = 0.5
        self.net['pore.area'] = 1.0
        self.net['throat.cross_sectional_area'] = 1.0
        self.phase = op.phase.Phase(network=self.net)
        self.phase['pore.diffusivity'] = 1.3
        self.phase['pore.molecular_weight'] = 0.029
        self.phase['pore.temperature'] = 345.0
        self.size_factors = np.ones([self.net.Nt, 3])*(0.123, 0.981, 0.551)

    def test_generic_diffusive_size_factors_as_dict(self):
        self.net['throat.diffusive_size_factors'] = self.size_factors
        mod = op.models.physics.diffusive_conductance.generic_diffusive
        self.phase.add_model(propname='throat.diffusive_conductance', model=mod)
        self.phase.regenerate_models()
        actual = self.phase['throat.diffusive_conductance'].mean()
        assert_allclose(actual, desired=0.091204832 * 1.3)
        del self.net["throat.diffusive_size_factors"]

    def test_generic_diffusive_size_factors_as_array(self):
        self.net['throat.diffusive_size_factors'] = 0.896
        mod = op.models.physics.diffusive_conductance.generic_diffusive
        self.phase.add_model(propname='throat.diffusive_conductance', model=mod)
        self.phase.regenerate_models()
        actual = self.phase['throat.diffusive_conductance'].mean()
        assert_allclose(actual, desired=0.896 * 1.3)
        del self.net["throat.diffusive_size_factors"]

    def test_generic_diffusive_partial_domain(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        phase = op.phase.Phase(network=net)
        phase.set_label('pore.left', pores=net['pore.left'])
        phase.set_label(label='left',
                        throats=net.find_neighbor_throats(net.pores('left')))
        phase['pore.diffusivity'] = 1.3
        net['throat.diffusive_size_factors'] = 0.5
        mod = op.models.physics.diffusive_conductance.generic_diffusive
        phase.add_model(propname='throat.diffusive_conductance',
                        model=mod,
                        domain='left',
                        regen_mode='deferred')
        phase.regenerate_models()
        actual = phase['throat.diffusive_conductance'].mean()
        # assert_allclose(actual, 0.65)
        # assert len(phase['throat.diffusive_conductance']) == 5

    def test_ordinary_diffusion_size_factors_given_as_dict(self):
        self.net['throat.diffusive_size_factors'] = self.size_factors
        mod = op.models.physics.diffusive_conductance.ordinary_diffusion
        self.phase.add_model(propname='throat.diffusive_conductance', model=mod)
        self.phase.regenerate_models()
        actual = self.phase['throat.diffusive_conductance'].mean()
        assert_allclose(actual, desired=0.091204832 * 1.3)
        del self.net["throat.diffusive_size_factors"]

    def test_ordinary_diffusion_size_factors_given_as_array(self):
        self.net['throat.diffusive_size_factors'] = 0.896
        mod = op.models.physics.diffusive_conductance.ordinary_diffusion
        self.phase.add_model(propname='throat.diffusive_conductance', model=mod)
        self.phase.regenerate_models()
        actual = self.phase['throat.diffusive_conductance'].mean()
        assert_allclose(actual, desired=0.896 * 1.3)
        del self.net["throat.diffusive_size_factors"]

    def test_mixed_diffusion(self):
        self.net['throat.diffusive_size_factors'] = self.size_factors
        mod = op.models.physics.diffusive_conductance.mixed_diffusion
        self.phase.add_model(propname='throat.diffusive_conductance', model=mod)
        self.phase.regenerate_models()
        actual = self.phase['throat.diffusive_conductance'].mean()
        assert_allclose(actual, desired=0.117568, rtol=1e-5)
        del self.net["throat.diffusive_size_factors"]

    # def test_multiphase_diffusion(self):
    #     net = op.network.Cubic(shape=[1, 6, 1])
    #     net['throat.diffusive_size_factors'] = self.size_factors_dict
    #     air = op.phase.Air(network=net)
    #     water = op.phase.Water(network=net)
    #     water['pore.diffusivity'] = 1e-11
    #     from openpnm import contrib
    #     m = contrib.MultiPhase(network=net, phases=[air, water])
    #     m.set_occupancy(phase=air, pores=[0, 1, 2])
    #     m.set_occupancy(phase=water, pores=[3, 4, 5])
    #     const = op.models.misc.constant
    #     K_water_air = 0.5
    #     m.set_binary_partition_coef(phases=[water, air], model=const, value=K_water_air)
    #     m._set_automatic_throat_occupancy()
    #     mdiff = op.models.physics.diffusive_conductance.multiphase_diffusion
    #     m.add_model(propname="throat.diffusive_conductance", model=mdiff)
    #     g = m["throat.diffusive_conductance"]
    #     # Diffusive conductance for MultiPhase must be Nt by 2 (not Nt by 1)
    #     assert g.shape == (net.Nt, 2)
    #     # Columns 1, 2 of conductance must be equal except for interface throat
    #     assert_allclose(g[:, 0][[0, 1, 3, 4]], g[:, 1][[0, 1, 3, 4]])
    #     # G12 and G21 at interface must differ (ratio must be K_water_air)
    #     assert_allclose(g[2, 0] / g[2, 1], 1 / K_water_air)
    #     assert_allclose(g.mean(), 7.64303956e-7, rtol=1e-5)

    def test_taylor_aris_diffusion(self):
        self.net['throat.diffusive_size_factors'] = self.size_factors
        self.phase['pore.pressure'] = np.linspace(0, 20, self.net.Np)
        self.phase['throat.hydraulic_conductance'] = 1
        mod = op.models.physics.diffusive_conductance.taylor_aris_diffusion
        self.phase.add_model(propname='throat.diffusive_conductance', model=mod)
        self.phase.regenerate_models()
        actual = np.array([
            self.phase['throat.diffusive_conductance'].mean(),
            self.phase['throat.diffusive_conductance'].max(),
            self.phase['throat.diffusive_conductance'].min()]
        )
        desired = np.array([0.121193, 0.126131, 0.118578])
        assert_allclose(actual, desired, rtol=1e-5)

    def test_spheres_and_cylinders(self):
        pn = op.network.Cubic(shape=[2, 1, 1])
        Dp1, Dp2 = 0.5, 0.25
        Dt = 0.1
        pn['pore.diameter'] = [Dp1, Dp2]
        pn['throat.diameter'] = Dt

        f = op.models.geometry.diffusive_size_factors.spheres_and_cylinders
        pn.add_model(propname='throat.diffusive_size_factors', model=f)
        Scalc = pn['throat.diffusive_size_factors'].flatten()
        # Pore 1
        Lp1 = (Dp1/2)*np.cos(np.arcsin(Dt/Dp1))
        Sp1 = 1/(np.pi*(Dp1/2))*np.arctanh(Lp1/(Dp1/2))
        assert_allclose(1/Sp1, Scalc[0], rtol=1e-7)
        # Pore 2
        Lp2 = (Dp2/2)*np.cos(np.arcsin(Dt/Dp2))
        Sp2 = 1/(np.pi*(Dp2/2))*np.arctanh(Lp2/(Dp2/2))
        assert_allclose(1/Sp2, Scalc[2], rtol=1e-7)
        # Throat
        Lt = 1 - Lp1 - Lp2
        St = Lt/(np.pi*(Dt/2)**2)
        assert_allclose(1/St, Scalc[1], rtol=1e-6)
        # Finally check conductance
        f = op.models.physics.diffusive_conductance.generic_diffusive
        phase = op.phase.Phase(network=pn)
        phase['pore.diffusivity'] = 1.0
        phase['throat.diffusivity'] = 1.0
        phase.add_model(propname='throat.diffusive_conductance', model=f)
        Gcalc = phase['throat.diffusive_conductance'][0]
        G = (Sp1 + St + Sp2)**-1
        assert_allclose(G, Gcalc, rtol=1e-7)

    def test_spheres_and_cylinders_oversized(self):
        pn = op.network.Cubic(shape=[2, 1, 1])
        Dp1, Dp2 = 1.75, 1.25
        Dt = 0.1
        pn['pore.diameter'] = [Dp1, Dp2]
        pn['throat.diameter'] = Dt

        f = op.models.geometry.diffusive_size_factors.spheres_and_cylinders
        pn.add_model(propname='throat.diffusive_size_factors', model=f)
        Scalc = pn['throat.diffusive_size_factors'].flatten()
        # Pore 1
        Lp1 = (Dp1/2)*np.cos(np.arcsin(Dt/Dp1))
        Sp1 = 1/(np.pi*(Dp1/2))*np.arctanh(Lp1/(Dp1/2))
        assert_allclose(1/Sp1, Scalc[0], rtol=1e-7)
        # Pore 2
        Lp2 = (Dp2/2)*np.cos(np.arcsin(Dt/Dp2))
        Sp2 = 1/(np.pi*(Dp2/2))*np.arctanh(Lp2/(Dp2/2))
        assert_allclose(1/Sp2, Scalc[2], rtol=1e-7)
        # Throat
        Lt = 1 - Lp1 - Lp2
        St = Lt/(np.pi*(Dt/2)**2)
        assert_allclose(1/St, Scalc[1], rtol=1e-6)
        # Finally check conductance
        f = op.models.physics.diffusive_conductance.generic_diffusive
        phase = op.phase.Phase(network=pn)
        phase['pore.diffusivity'] = 1.0
        phase['throat.diffusivity'] = 1.0
        phase.add_model(propname='throat.diffusive_conductance', model=f)
        Gcalc = phase['throat.diffusive_conductance'][0]
        G = (Sp1 + St + Sp2)**-1
        assert_allclose(G, Gcalc, rtol=1e-7)

    def test_cones_and_cylinders(self):
        pn = op.network.Cubic(shape=[2, 1, 1])
        Dp1, Dp2 = 0.5, 0.25
        Dt = 0.1
        pn['pore.diameter'] = [Dp1, Dp2]
        pn['throat.diameter'] = Dt

        f = op.models.geometry.diffusive_size_factors.cones_and_cylinders
        pn.add_model(propname='throat.diffusive_size_factors', model=f)
        Scalc = pn['throat.diffusive_size_factors'].flatten()
        # Pore 1
        Sp1 = (Dp1/2)/(np.pi*(Dp1/2*Dt/2))
        assert_allclose(1/Sp1, Scalc[0], rtol=1e-7)
        # Pore 2
        Sp2 = (Dp2/2)/(np.pi*(Dp2/2*Dt/2))
        assert_allclose(1/Sp2, Scalc[2], rtol=1e-7)
        # Throat
        Lt = 1 - Dp1/2 - Dp2/2
        St = Lt/(np.pi*(Dt/2)**2)
        assert_allclose(1/St, Scalc[1], rtol=1e-6)
        # Finally check conductance
        f = op.models.physics.diffusive_conductance.generic_diffusive
        phase = op.phase.Phase(network=pn)
        phase['pore.diffusivity'] = 1.0
        phase['throat.diffusivity'] = 1.0
        phase.add_model(propname='throat.diffusive_conductance', model=f)
        Gcalc = phase['throat.diffusive_conductance'][0]
        G = (Sp1 + St + Sp2)**-1
        assert_allclose(G, Gcalc, rtol=1e-7)

    def test_intersecting_cones(self):
        pn = op.network.Cubic(shape=[2, 1, 1])
        Dp1, Dp2 = 0.5, 0.25
        Dt = 0.1
        pn['pore.diameter'] = [Dp1, Dp2]
        pn['throat.diameter'] = Dt
        pn['throat.coords'] = [[1.0, 0.5, 0.5]]

        f = op.models.geometry.diffusive_size_factors.intersecting_cones
        pn.add_model(propname='throat.diffusive_size_factors', model=f)
        Scalc = pn['throat.diffusive_size_factors'].flatten()
        # Pore 1
        Lp1 = 0.5
        Sp1 = (Lp1)/(np.pi*(Dp1/2*Dt/2))
        assert_allclose(1/Sp1, Scalc[0], rtol=1e-7)
        # Pore 2
        Lp2 = 0.5
        Sp2 = (Lp2)/(np.pi*(Dp2/2*Dt/2))
        assert_allclose(1/Sp2, Scalc[2], rtol=1e-7)
        # Throat
        St = 0
        # Finally check conductance
        f = op.models.physics.diffusive_conductance.generic_diffusive
        phase = op.phase.Phase(network=pn)
        phase['pore.diffusivity'] = 1.0
        phase['throat.diffusivity'] = 1.0
        phase.add_model(propname='throat.diffusive_conductance', model=f)
        Gcalc = phase['throat.diffusive_conductance'][0]
        G = (Sp1 + St + Sp2)**-1
        assert_allclose(G, Gcalc, rtol=1e-7)













if __name__ == '__main__':

    t = DiffusiveConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
