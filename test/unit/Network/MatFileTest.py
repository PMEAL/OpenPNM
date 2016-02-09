import OpenPNM
import scipy as sp
import os


class MatFileTest:
    def setup_class(self):
        pass

    def test_loading(self):
        fname = os.path.join(FIXTURE_DIR, 'example_network.mat')
        pn = OpenPNM.Network.MatFile(filename=fname)
        assert pn.Np == 105
        assert pn.Nt == 238
        assert pn.props() == ['throat.conns', 'pore.coords']

    def test_extra_args(self):
        fname = os.path.join(FIXTURE_DIR, 'example_network.mat')
        pn = OpenPNM.Network.MatFile(filename=fname,
                                     xtra_pore_data='type',
                                     xtra_throat_data='type')
        assert sp.all(sp.unique(pn['pore.type']) == [0, 1, 6])
