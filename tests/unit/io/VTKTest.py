import py
import os
import numpy as np
import openpnm as op
from pathlib import Path


class VTKTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True
        self.net = op.network.Cubic(shape=[2, 2, 2])
        self.net['pore.boo'] = 1
        self.net['throat.boo'] = 1

        self.phase_1 = op.phase.GenericPhase(network=self.net)
        self.phase_1['pore.bar'] = 2
        self.phase_1['throat.bar'] = 2
        self.phase_2 = op.phase.GenericPhase(network=self.net)
        self.phase_2['pore.bar'] = 2
        self.phase_2['throat.bar'] = 2

        self.net['pore.object'] = np.ones(self.net.Np, dtype=object)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_save_network(self, tmpdir):
        fname = Path(tmpdir,  'test_save_vtk_1.vtp')
        op.io.to_vtk(network=self.net, filename=fname)
        assert fname.is_file()
        os.remove(fname)

    def test_save_network_and_phases(self, tmpdir):
        fname = Path(tmpdir,  'test_save_vtk_2.vtp')
        op.io.to_vtk(network=self.net, phases=self.phase_1, filename=fname)
        assert fname.is_file()
        os.remove(fname)


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = VTKTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
