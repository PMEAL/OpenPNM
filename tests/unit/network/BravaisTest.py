import openpnm as op
import pytest
ws = op.Workspace()


class BravaisTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        ws.clear()

    def test_generation_fcc(self):
        fcc = op.network.Bravais(shape=[2, 2, 2], mode='fcc')
        assert fcc.Np == 63
        assert fcc.Nt == 294
        assert fcc.num_pores('pore.corner_sites') == 27
        assert fcc.num_pores('pore.face_sites') == 36
        assert fcc.num_throats('throat.corner_to_corner') == 54
        assert fcc.num_throats('throat.corner_to_face') == 240

    def test_fcc_add_boundary_pores(self):
        fcc = op.network.Bravais(shape=[2, 2, 2], mode='fcc')
        Np = fcc.Np
        op.topotools.add_boundary_pores(network=fcc, pores=fcc.pores('left'),
                                        apply_label='left_boundary')
        op.topotools.add_boundary_pores(network=fcc, pores=fcc.pores('right'),
                                        apply_label='right_boundary')
        assert fcc.Np > Np
        assert fcc.Np == (Np + 18 + 8)
        assert 'pore.left_boundary' in fcc.keys()
        assert 'pore.right_boundary' in fcc.keys()
        assert 'pore.top_boundary' not in fcc.keys()
        assert 'pore.bottom_boundary' not in fcc.keys()

    def test_generation_bcc(self):
        bcc = op.network.Bravais(shape=[2, 2, 2], mode='bcc')
        assert bcc.Np == 35
        assert bcc.Nt == 130
        assert bcc.num_pores('pore.corner_sites') == 27
        assert bcc.num_pores('pore.body_sites') == 8
        assert bcc.num_throats('throat.corner_to_corner') == 54
        assert bcc.num_throats('throat.corner_to_body') == 64
        assert bcc.num_throats('throat.body_to_body') == 12

    def test_bcc_add_boundary_pores(self):
        bcc = op.network.Bravais(shape=[2, 2, 2], mode='bcc')
        Np = bcc.Np
        op.topotools.add_boundary_pores(network=bcc, pores=bcc.pores('left'),
                                        apply_label='left_boundary')
        op.topotools.add_boundary_pores(network=bcc, pores=bcc.pores('right'),
                                        apply_label='right_boundary')
        assert bcc.Np > Np
        assert bcc.Np == (Np + 18)
        assert 'pore.left_boundary' in bcc.keys()
        assert 'pore.right_boundary' in bcc.keys()
        assert 'pore.top_boundary' not in bcc.keys()
        assert 'pore.bottom_boundary' not in bcc.keys()

    def test_generation_hcp(self):
        with pytest.raises(NotImplementedError):
            op.network.Bravais(shape=[3, 3, 3], mode='hcp')

    def test_generation_bad_mode(self):
        with pytest.raises(Exception):
            op.network.Bravais(shape=[3, 3, 3], mode='bad_mode')

    def test_generation_sc(self):
        sc = op.network.Bravais(shape=[3, 3, 3], mode='sc')
        assert sc.Np == 27
        assert sc.Nt == 54

    def test_generation_too_small(self):
        with pytest.raises(Exception):
            op.network.Bravais(shape=[1, 2, 3], mode='')


if __name__ == '__main__':

    t = BravaisTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
