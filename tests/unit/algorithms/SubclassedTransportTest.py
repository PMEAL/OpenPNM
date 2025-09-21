from numpy.testing import assert_allclose

import openpnm as op


class SubclassedTransportTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[9, 9, 9])
        self.net.add_model_collection(
            op.models.collections.geometry.spheres_and_cylinders)
        self.net.regenerate_models()
        self.phase = op.phase.Phase(network=self.net)
        self.phase.add_model_collection(
            op.models.collections.physics.standard)
        self.phase.regenerate_models()

    def test_fickian_diffusion(self):
        alg = op.algorithms.FickianDiffusion(network=self.net,
                                             phase=self.phase)
        self.phase['throat.diffusive_conductance'] = 1
        Pin = self.net.pores('top')
        Pout = self.net.pores('bottom')
        alg.set_value_BC(pores=Pin, values=1)
        alg.set_value_BC(pores=Pout, values=0)
        alg.run()
        rate_in = alg.rate(pores=Pin)[0]
        rate_out = alg.rate(pores=Pout)[0]
        rate_total = alg.rate(self.net.pores("surface"))[0]
        assert_allclose(rate_in, 10.125)
        assert_allclose(rate_in, -rate_out)
        assert rate_in > 0
        assert abs(rate_total) < rate_in * 1e-13

    def test_stokes_flow(self):
        alg = op.algorithms.StokesFlow(network=self.net, phase=self.phase)
        self.phase['throat.hydraulic_conductance'] = 1
        Pin = self.net.pores('top')
        Pout = self.net.pores('bottom')
        alg.set_value_BC(pores=Pin, values=101325)
        alg.set_value_BC(pores=Pout, values=0)
        alg.run()
        rate_in = alg.rate(pores=Pin)[0]
        rate_out = alg.rate(pores=Pout)[0]
        rate_total = alg.rate(self.net.pores("surface"))[0]
        assert_allclose(rate_in, 1025915.625)
        assert_allclose(rate_in, -rate_out)
        assert rate_in > 0
        assert abs(rate_total) < rate_in * 1e-13

    def test_forurier_conduction(self):
        alg = op.algorithms.FourierConduction(network=self.net,
                                              phase=self.phase)
        self.phase['throat.thermal_conductance'] = 1
        Pin = self.net.pores('top')
        Pout = self.net.pores('bottom')
        alg.set_value_BC(pores=Pin, values=273)
        alg.set_value_BC(pores=Pout, values=0)
        alg.run()
        rate_in = alg.rate(pores=Pin)[0]
        rate_out = alg.rate(pores=Pout)[0]
        rate_total = alg.rate(self.net.pores("surface"))[0]
        assert_allclose(rate_in, 2764.125)
        assert_allclose(rate_in, -rate_out)
        assert rate_in > 0
        assert abs(rate_total) < rate_in * 1e-13

    def test_ohmic_conduction(self):
        alg = op.algorithms.OhmicConduction(network=self.net, phase=self.phase)
        self.phase['throat.electrical_conductance'] = 1
        Pin = self.net.pores('top')
        Pout = self.net.pores('bottom')
        alg.set_value_BC(pores=Pin, values=1)
        alg.set_value_BC(pores=Pout, values=0)
        alg.run()
        rate_in = alg.rate(pores=Pin)[0]
        rate_out = alg.rate(pores=Pout)[0]
        rate_total = alg.rate(self.net.pores("surface"))[0]
        assert_allclose(rate_in, 10.125)
        assert_allclose(rate_in, -rate_out)
        assert rate_in > 0
        assert abs(rate_total) < rate_in * 1e-13

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = SubclassedTransportTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()
