import py
import os
import openpnm as op
from pathlib import Path


class MARockTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.clear()

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_load_MARock(self, tmpdir):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/3DMA-Castlegate')
        project = op.io.from_marock(path=path)
        assert len(project) == 1
        net = project.network
        assert net.Np == 9915
        assert net.Nt == 21805
        a = {'pore.ID_number', 'pore.boundary_type', 'pore.coordination',
             'pore.coords', 'pore.volume', 'throat.area', 'throat.conns',
             'throat.coords'}
        assert a.issubset(net.props())


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = MARockTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
