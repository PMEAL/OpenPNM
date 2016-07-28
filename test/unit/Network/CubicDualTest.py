import OpenPNM as op
import scipy as sp


class CubicDualTest:

    def test_generation(self):
        net = op.Network.CubicDual(shape=[5, 5, 5], label_1='primary',
                                   label_2='secondary')
        assert net.Np == 285
        assert net.Nt == 1052
        assert net.num_pores('primary') == 125
        assert net.num_pores('secondary') == 160
        assert net.num_pores('surface') == 250
        assert net.num_throats('primary') == 300
        assert net.num_throats('secondary') == 240
        assert net.num_throats('surface') == 96
        assert net.num_throats('interconnect') == 512

    def test_add_boundary_pores(self):
        net = op.Network.CubicDual(shape=[5, 5, 5], label_1='primary',
                                   label_2='secondary')
        Ps = net.pores(labels=['surface', 'bottom'], mode='intersection')
        net.add_boundary_pores(pores=Ps, offset=[0, 0, -0.5])
        Ps2 = net.pores(labels=['boundary'], mode='intersection')
        assert Ps.size == Ps2.size
        assert ~sp.any(sp.in1d(Ps, Ps2))
