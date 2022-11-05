import os
import numpy as np
import scipy as sp
import openpnm as op
import networkx as nx
from pathlib import Path


class StatoilTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_load_F42A(self, tmpdir):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/ICL-SandPack(F42A)')
        net = op.io.network_from_statoil(path=path.resolve(), prefix='F42A')
        self.net = net
        assert net.Np == 1246
        assert net.Nt == 2654
        assert np.shape(net['pore.coords']) == (1246, 3)
        assert np.shape(net['throat.conns']) == (2654, 2)
        assert 'pore.radius' in net.keys()

    def test_load_Berea(self, tmpdir):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/ICL-Sandstone(Berea)')
        net = op.io.network_from_statoil(path=path, prefix='Berea')
        assert net.Np == 6298
        assert net.Nt == 12098
        assert np.shape(net['pore.coords']) == (6298, 3)
        assert np.shape(net['throat.conns']) == (12098, 2)
        assert 'pore.radius' in net.keys()
        assert np.all(net.find_neighbor_pores(pores=1000) == [221, 1214])


if __name__ == '__main__':
    import py
    # All the tests in this file can be run with 'playing' this file
    t = StatoilTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('Running test: {item}')
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
