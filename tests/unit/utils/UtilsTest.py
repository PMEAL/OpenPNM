import openpnm as op
import scipy as sp
import pytest


class UtilsTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_tic_toc(self):
        with pytest.raises(Exception):
            op.utils.toc()
        op.utils.tic()
        t1 = op.utils.toc()
        assert t1 is None
        t2 = op.utils.toc(quiet=True)
        assert t2 >= 0

    def test_nested_dict(self):
        d = op.utils.NestedDict()
        d['top']['middle']['bottom'] = 1
        assert d == {'top': {'middle': {'bottom': 1}}}
        s = d.__str__()
        assert s == '-top\n--middle\n---bottom\n'
        a = d.to_dict()
        assert type(a) is dict

    def test_printable_list(self):
        L = op.utils.PrintableList(['item1', 'item2', 'item2'])
        s = L.__str__().split('\n')
        assert len(s) == 5

    def test_printable_dict(self):
        D = op.utils.PrintableDict(**{'item1': 1, 'item2': 2,
                                      'item3': sp.array([1, 2])})
        s = D.__str__().split('\n')
        assert len(s) == 7
        r = D.__repr__()
        assert r == "{'item1': 1, 'item2': 2, 'item3': array([1, 2])}"

    def test_is_symmetric_w_rtol(self):
        A = sp.array([[1, 2, 3], [2, 4, 6], [3.000001, 6, 99]])
        is_sym = op.utils.misc.is_symmetric(A, rtol=1e-4)
        assert is_sym
        is_sym = op.utils.misc.is_symmetric(A, rtol=1e-6)
        assert not is_sym

    def test_is_symmetric_FickianDiffusion_must_be_symmetric(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        geom = op.geometry.StickAndBall(network=net)
        air = op.phases.Air(network=net)
        _ = op.physics.Standard(network=net, phase=air, geometry=geom)
        fd = op.algorithms.FickianDiffusion(network=net, phase=air)
        fd.set_value_BC(pores=net.pores("left"), values=1.0)
        fd.set_value_BC(pores=net.pores("right"), values=0.1)
        assert op.utils.misc.is_symmetric(fd.A)

    def test_is_symmetric_AdvectionDiffusion_must_be_nonsymmetric(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        geom = op.geometry.StickAndBall(network=net)
        air = op.phases.Air(network=net)
        phys = op.physics.Standard(network=net, phase=air, geometry=geom)
        ad = op.algorithms.AdvectionDiffusion(network=net, phase=air)
        ad.set_value_BC(pores=net.pores("left"), values=1.0)
        ad.set_value_BC(pores=net.pores("right"), values=0.1)
        # Uniform pressure field --> no advection --> still symmetric
        assert op.utils.misc.is_symmetric(ad.A)
        sf = op.algorithms.StokesFlow(network=net, phase=air)
        sf.set_value_BC(pores=net.pores("left"), values=1.0)
        sf.set_value_BC(pores=net.pores("right"), values=0.0)
        sf.run()
        air.update(sf.results())
        phys.regenerate_models()
        ad.settings.update({"cache_A": False, "cache_b": False})
        ad._build_A()
        # Non-uniform pressure field --> positive advection --> non-symmetric
        assert not op.utils.misc.is_symmetric(ad.A)


if __name__ == '__main__':

    t = UtilsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
