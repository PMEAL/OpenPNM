import openpnm as op
ws = op.Workspace()


class LabelTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        ws.clear()

    def test_bravais_fcc(self):
        net = op.network.Bravais(shape=[4, 4, 4], mode='fcc')
        assert 'pore.surface' in net.keys()
        assert 'pore.internal' in net.keys()
        assert net.num_pores('left') == 25
        assert net.num_pores(['left', 'top'], mode='xnor') == 4

    def test_bravais_bcc(self):
        net = op.network.Bravais(shape=[4, 4, 4], mode='bcc')
        assert 'pore.surface' in net.keys()
        assert 'pore.internal' in net.keys()
        assert net.num_pores('left') == 16
        assert net.num_pores(['left', 'top'], mode='xnor') == 4

    def test_bravais_sc(self):
        net = op.network.Bravais(shape=[4, 4, 4], mode='sc')
        assert 'pore.surface' in net.keys()
        assert 'pore.internal' in net.keys()
        assert net.num_pores('left') == 16
        assert net.num_pores(['left', 'top'], mode='xnor') == 4

    def test_cubic_2D(self):
        net = op.network.Cubic(shape=[5, 5, 1])
        assert 'pore.top' not in net.labels()
        assert 'pore.bottom' not in net.labels()
        assert net.num_pores('surface') == 16
        net = op.network.Cubic(shape=[5, 1, 5])
        assert 'pore.left' not in net.labels()
        assert 'pore.right' not in net.labels()
        assert net.num_pores('surface') == 16
        net = op.network.Cubic(shape=[1, 5, 5])
        assert 'pore.front' not in net.labels()
        assert 'pore.back' not in net.labels()
        assert net.num_pores('surface') == 16

if __name__ == '__main__':

    t = LabelTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
