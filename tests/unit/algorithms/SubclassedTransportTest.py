import openpnm as op
from numpy.testing import assert_approx_equal


class SubclassedTransportTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[9, 9, 9])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phase['pore.viscosity'] = 1e-3

    def test_fickian_diffusion(self):
        alg = op.algorithms.FickianDiffusion(network=self.net,
                                             phase=self.phase)
        self.phys['throat.diffusive_conductance'] = 1
        Pin = self.net.pores('top')
        Pout = self.net.pores('bottom')
        alg.set_value_BC(pores=Pin, values=1)
        alg.set_value_BC(pores=Pout, values=0)
        alg.run()
        Deff = alg.calc_effective_diffusivity()
        Deff = alg.calc_effective_diffusivity(domain_area=81, domain_length=9)
        Deff = alg.calc_effective_diffusivity(domain_area=81, domain_length=9,
                                              inlets=Pin, outlets=Pout)
        assert_approx_equal(Deff, 1.12500)

    def test_stokes_flow(self):
        alg = op.algorithms.StokesFlow(network=self.net, phase=self.phase)
        self.phys['throat.hydraulic_conductance'] = 1
        Pin = self.net.pores('top')
        Pout = self.net.pores('bottom')
        alg.set_value_BC(pores=Pin, values=101325)
        alg.set_value_BC(pores=Pout, values=0)
        alg.run()
        Keff = alg.calc_effective_permeability()
        Keff = alg.calc_effective_permeability(domain_area=81, domain_length=9)
        Keff = alg.calc_effective_permeability(domain_area=81, domain_length=9,
                                               inlets=Pin, outlets=Pout)
        assert_approx_equal(Keff, 0.001125)

    def test_forurier_conduction(self):
        alg = op.algorithms.FourierConduction(network=self.net,
                                              phase=self.phase)
        self.phys['throat.thermal_conductance'] = 1
        Pin = self.net.pores('top')
        Pout = self.net.pores('bottom')
        alg.set_value_BC(pores=Pin, values=101325)
        alg.set_value_BC(pores=Pout, values=0)
        alg.run()
        Keff = alg.calc_effective_conductivity()
        Keff = alg.calc_effective_conductivity(domain_area=81, domain_length=9)
        Keff = alg.calc_effective_conductivity(domain_area=81, domain_length=9,
                                               inlets=Pin, outlets=Pout)
        assert_approx_equal(Keff, 1.125)

    def test_ohmic_conduction(self):
        alg = op.algorithms.OhmicConduction(network=self.net, phase=self.phase)
        self.phys['throat.electrical_conductance'] = 1
        Pin = self.net.pores('top')
        Pout = self.net.pores('bottom')
        alg.set_value_BC(pores=Pin, values=101325)
        alg.set_value_BC(pores=Pout, values=0)
        alg.run()
        Keff = alg.calc_effective_conductivity()
        Keff = alg.calc_effective_conductivity(domain_area=81, domain_length=9)
        Keff = alg.calc_effective_conductivity(domain_area=81, domain_length=9,
                                               inlets=Pin, outlets=Pout)
        assert_approx_equal(Keff, 1.125)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = SubclassedTransportTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
