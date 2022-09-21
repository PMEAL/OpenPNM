import openpnm as op
import pytest
ws = op.Workspace()


class BravaisTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        ws.clear()

    def test_generation_fcc(self):
        fcc = op.network.FaceCenteredCubic(shape=[3, 3, 3])
        assert fcc.Np == 63
        assert fcc.Nt == 294
        assert fcc.num_pores('pore.corner') == 27
        assert fcc.num_pores('pore.face') == 36
        assert fcc.num_throats('throat.corner_to_corner') == 54
        assert fcc.num_throats('throat.corner_to_face') == 240

    def test_fcc_add_boundary_pores(self):
        fcc = op.network.FaceCenteredCubic(shape=[3, 3, 3])
        Np = fcc.Np
        fcc.add_boundary_pores(labels=['left', 'right'], spacing=1)
        assert fcc.Np > Np
        assert fcc.Np == (Np + 18 + 8)
        assert 'pore.left_boundary' in fcc.keys()
        assert 'pore.right_boundary' in fcc.keys()
        assert 'pore.top_boundary' not in fcc.keys()
        assert 'pore.bottom_boundary' not in fcc.keys()

    def test_generation_bcc(self):
        bcc = op.network.BodyCenteredCubic(shape=[3, 3, 3])
        assert bcc.Np == 35
        assert bcc.Nt == 130
        assert bcc.num_pores('pore.corner') == 27
        assert bcc.num_pores('pore.body') == 8
        assert bcc.num_throats('throat.corner_to_corner') == 54
        assert bcc.num_throats('throat.corner_to_body') == 64
        assert bcc.num_throats('throat.body_to_body') == 12

    def test_bcc_add_boundary_pores(self):
        bcc = op.network.BodyCenteredCubic(shape=[3, 3, 3])
        Np = bcc.Np
        bcc.add_boundary_pores(labels=['left', 'right'], spacing=1)
        assert bcc.Np > Np
        assert bcc.Np == (Np + 18)
        assert 'pore.left_boundary' in bcc.keys()
        assert 'pore.right_boundary' in bcc.keys()
        assert 'pore.top_boundary' not in bcc.keys()
        assert 'pore.bottom_boundary' not in bcc.keys()

    def test_generation_too_small(self):
        with pytest.raises(Exception):
            op.network.BodyCentered(shape=[1, 2, 3], mode='')


if __name__ == '__main__':

    t = BravaisTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
