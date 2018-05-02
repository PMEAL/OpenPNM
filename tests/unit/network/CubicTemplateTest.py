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


if __name__ == '__main__':

    t = CubicTemplateTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
