import os
import numpy as np
import openpnm as op
from pathlib import Path


class iMorphTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_load(self):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/iMorph-Sandstone')
        project = op.io.iMorph.import_data(path)
        assert len(project) == 1
        net = project.network
        assert net.Np == 1518
        assert net.Nt == 2424
        assert np.shape(net['pore.coords']) == (1518, 3)
        assert np.shape(net['throat.conns']) == (2424, 2)
        a = {'pore.volume', 'pore.types', 'throat.volume', 'throat.types'}
        assert a.issubset(net.props())
        a = {'pore.internal', 'pore.top_boundary', 'pore.bottom_boundary',
             'pore.front_boundary', 'pore.back_boundary', 'pore.left_boundary',
             'pore.right_boundary'}
        assert a.issubset(net.labels())


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = iMorphTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
