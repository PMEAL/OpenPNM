import numpy as np
import openpnm as op
from skimage.morphology import ball, disk


class CubicTemplateTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_2D_template(self):
        net = op.network.CubicTemplate(template=disk(10), spacing=1)
        assert net.Np == 317
        assert net.Nt == 592

    def test_3D_template(self):
        net = op.network.CubicTemplate(template=ball(5), spacing=1)
        assert net.Np == 515
        assert net.Nt == 1302

    # def test_labels(self):
    #     template = np.array(
    #         [[1, 1, 1, 1, 1],
    #          [1, 1, 0, 1, 1],
    #          [1, 1, 0, 0, 1],
    #          [1, 0, 0, 0, 1],
    #          [1, 1, 0, 1, 1]]
    #     )
    #     net = op.network.CubicTemplate(template=template)
    #     # Test "surface" label
    #     Ps_surf_desired = np.array([0, 1, 2, 3, 4, 5, 8, 9, 11, 12, 13, 14, 15, 16, 17])
    #     Ps_surf = net.pores("surface")
    #     np.testing.assert_allclose(Ps_surf, Ps_surf_desired)
    #     # Test "internal_surface" label
    #     Ps_int_surf_desired = np.array([6, 7, 10])
    #     Ps_int_surf = net.pores("internal_surface")
    #     np.testing.assert_allclose(Ps_int_surf, Ps_int_surf_desired)


if __name__ == '__main__':

    t = CubicTemplateTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
