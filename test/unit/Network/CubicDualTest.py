import OpenPNM as op


class CubicDualTest:

    def test_generation(self):
        net = op.Network.CubicDual(shape=[5, 5, 5], label_1='primary',
                                   label_2='secondary')
        assert net.Np == 285
        assert net.Nt == 1052
        assert net.num_pores('primary') == 125
        assert net.num_pores('secondary') == 125
        assert net.num_pores('surface') == 250
        assert net.num_throats('primary') == 300
        assert net.num_throats('secondary') == 240
        assert net.num_throats('surface') == 96
        assert net.num_throats('interconnect') == 512
