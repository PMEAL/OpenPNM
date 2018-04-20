import openpnm as op
import scipy as sp
import pytest


class SubclassedTransportTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[9, 9, 9])
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)
        self.phase = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)

    def test_fickian_diffusion(self):
        alg = op.algorithms.FickianDiffusion(network=self.net,
                                             phase=self.phase)
        self.phys['throat.diffusive_conductance'] = 1
        alg.set_dirichlet_BC(pores=self.net.pores('top'), values=1)
        alg.set_dirichlet_BC(pores=self.net.pores('bottom'), values=0)
        alg.run()
        alg.domain_area = 81
        alg.domain_length = 9
        Deff = alg.calc_eff_diffusivity()
        assert sp.around(Deff, decimals=10) == 0.0275097564

    def test_stokes_flow(self):
        alg = op.algorithms.StokesFlow(network=self.net, phase=self.phase)
        self.phys['throat.hydraulic_conductance'] = 1
        alg.set_dirichlet_BC(pores=self.net.pores('top'), values=101325)
        alg.set_dirichlet_BC(pores=self.net.pores('bottom'), values=0)
        alg.run()
        alg.domain_area = 81
        alg.domain_length = 9
        Keff = alg.calc_eff_permeability()
        assert sp.around(Keff, decimals=10) == 2.075e-05

    def test_forurier_conduction(self):
        alg = op.algorithms.FourierConduction(network=self.net,
                                              phase=self.phase)
        self.phys['throat.thermal_conductance'] = 1
        alg.set_dirichlet_BC(pores=self.net.pores('top'), values=101325)
        alg.set_dirichlet_BC(pores=self.net.pores('bottom'), values=0)
        alg.run()
        alg.domain_area = 81
        alg.domain_length = 9
        Keff = alg.calc_effective_conductivity()
        assert sp.around(Keff, decimals=10) == 1.125

    def test_ohmic_conduction(self):
        alg = op.algorithms.OhmicConduction(network=self.net, phase=self.phase)
        self.phys['throat.electrical_conductance'] = 1
        alg.set_dirichlet_BC(pores=self.net.pores('top'), values=101325)
        alg.set_dirichlet_BC(pores=self.net.pores('bottom'), values=0)
        alg.run()
        alg.domain_area = 81
        alg.domain_length = 9
        Keff = alg.calc_effective_conductivity()
        assert sp.around(Keff, decimals=10) == 1.125

    def teardown_class(self):
        ws = op.core.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = SubclassedTransportTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
