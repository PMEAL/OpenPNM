import os
from pathlib import Path

import openpnm as op


class MARockTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.clear()

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_load_MARock(self):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/3DMA-Castlegate')
        net = op.io.network_from_marock(filename=path)
        assert hasattr(net, 'conns')
        assert net.Np == 9915
        assert net.Nt == 21805
        a = {'pore.ID_number', 'pore.boundary_type', 'pore.coordination',
             'pore.coords', 'pore.volume', 'throat.conns',
             'throat.coords', 'throat.cross_sectional_area'}
        assert a.issubset(net.props())


if __name__ == '__main__':
    import py

    # All the tests in this file can be run with 'playing' this file
    t = MARockTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()
