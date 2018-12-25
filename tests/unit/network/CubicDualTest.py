import openpnm as op


class CubicDualTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_generation_3D(self):
        net = op.network.CubicDual(shape=[5, 5, 5], label_1='primary',
                                   label_2='secondary')
        assert net.Np == 285
        assert net.Nt == 1436
        assert net.num_pores('all') == 285
        assert net.num_pores('back') == 41
        assert net.num_pores('bottom') == 41
        assert net.num_pores('front') == 41
        assert net.num_pores('internal') == 285
        assert net.num_pores('left') == 41
        assert net.num_pores('primary') == 125
        assert net.num_pores('right') == 41
        assert net.num_pores('secondary') == 160
        assert net.num_pores('surface') == 194
        assert net.num_pores('top') == 41
        assert net.num_throats('all') == 1436
        assert net.num_throats('interconnect') == 896
        assert net.num_throats('internal') == 1436
        assert net.num_throats('primary') == 300
        assert net.num_throats('secondary') == 240
        assert net.num_throats('surface') == 576
        assert net.num_throats('top') == 104
        assert net.num_throats('bottom') == 104
        assert net.num_throats('left') == 104
        assert net.num_throats('right') == 104
        assert net.num_throats('front') == 104
        assert net.num_throats('back') == 104

    def test_generation_2D_XY(self):
        net = op.network.CubicDual(shape=[5, 5, 1], label_1='primary',
                                   label_2='secondary')
        assert net.Np == 57
        assert net.Nt == 176
        assert net.num_pores('left') == 9
        assert net.num_pores('right') == 9
        assert net.num_pores('front') == 9
        assert net.num_pores('back') == 9
        assert net.num_throats('interconnect') == 96

    def test_generation_2D_XZ(self):
        net = op.network.CubicDual(shape=[5, 1, 5], label_1='primary',
                                   label_2='secondary')
        assert net.Np == 57
        assert net.Nt == 176
        assert net.num_pores('front') == 9
        assert net.num_pores('back') == 9
        assert net.num_pores('top') == 9
        assert net.num_pores('bottom') == 9
        assert net.num_throats('interconnect') == 96

    def test_generation_2D_YZ(self):
        net = op.network.CubicDual(shape=[1, 5, 5], label_1='primary',
                                   label_2='secondary')
        assert net.Np == 57
        assert net.Nt == 176
        assert net.num_pores('left') == 9
        assert net.num_pores('right') == 9
        assert net.num_pores('top') == 9
        assert net.num_pores('bottom') == 9
        assert net.num_throats('interconnect') == 96

    def test_generation_2D_2_dims(self):
        net = op.network.CubicDual(shape=[5, 5], label_1='primary',
                                   label_2='secondary')
        assert net.Np == 57
        assert net.Nt == 176
        assert net.num_pores('left') == 9
        assert net.num_pores('right') == 9
        assert net.num_pores('front') == 9
        assert net.num_pores('back') == 9
        assert net.num_throats('interconnect') == 96


if __name__ == '__main__':

    t = CubicDualTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
