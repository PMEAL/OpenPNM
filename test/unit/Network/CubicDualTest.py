import OpenPNM as op
import scipy as sp


class CubicDualTest:

    def test_generation(self):
        net = op.Network.CubicDual(shape=[5, 5, 5], label_1='primary',
                                   label_2='secondary')
        assert net.Np == 285
        assert net.Nt == 1436
        assert net.num_pores('all') == 491
        assert net.num_pores('back') == 61
        assert net.num_pores('bottom') == 61
        assert net.num_pores('front') == 61
        assert net.num_pores('internal') == 189
        assert net.num_pores('left') == 61
        assert net.num_pores('primary') == 275
        assert net.num_pores('right') == 61
        assert net.num_pores('secondary') == 216
        assert net.num_pores('surface') == 302
        assert net.num_pores('top') == 61
        assert net.num_throats('all') == 2590
        assert net.num_throats('interconnect') == 1600
        assert net.num_throats('internal') == 1690
        assert net.num_throats('primary') == 450
        assert net.num_throats('secondary') == 540
        assert net.num_throats('surface') == 900

    def test_add_boundary_pores(self):
        net = op.Network.CubicDual(shape=[5, 5, 5], label_1='primary',
                                   label_2='secondary')
        Ps = net.pores(labels=['surface', 'bottom'], mode='intersection')
        net.add_boundary_pores(pores=Ps, offset=[0, 0, -0.5])
        Ps2 = net.pores(labels=['boundary'], mode='intersection')
        assert Ps.size == Ps2.size
        assert ~sp.any(sp.in1d(Ps, Ps2))
