import os

import numpy as np

import openpnm as op


class XDMFTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True
        self.net = op.network.Cubic(shape=[2, 2, 2])
        self.net['pore.boo'] = 1
        self.net['throat.boo'] = 1

        self.phase_1 = op.phase.Phase(network=self.net)
        self.phase_1['pore.bar'] = 2
        self.phase_1['throat.bar'] = 2
        self.phase_2 = op.phase.Phase(network=self.net)
        self.phase_2['pore.bar'] = 2
        self.phase_2['throat.bar'] = 2
        # The follow is confirmation that bug 1456 is fixed
        self.phase_2['throat.blah.de.blah'] = 1
        self.phase_2['throat.blah.bloo'] = 1

        self.net['pore.object'] = np.ones(self.net.Np, dtype=object)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_save(self, tmpdir):
        fname = tmpdir.join('test_file')
        op.io.project_to_xdmf(self.net.project, filename=fname)
        os.remove(tmpdir.join('test_file.hdf'))
        os.remove(tmpdir.join('test_file.xmf'))


if __name__ == '__main__':
    import py

    # All the tests in this file can be run with 'playing' this file
    t = XDMFTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
