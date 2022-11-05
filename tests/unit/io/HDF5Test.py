import os
import numpy as np
import openpnm as op


class HDF5Test:

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

        self.net['pore.object'] = np.ones(self.net.Np, dtype=object)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_project_to_hdf5(self, tmpdir):
        fname = tmpdir.join(self.net.project.name)
        f = op.io.project_to_hdf5(project=self.net.project, filename=fname)
        assert list(f.keys()) == [self.net.name, self.phase_1.name, self.phase_2.name]
        filename = f.filename
        f.close()
        os.remove(filename)

    def test_print_hdf5(self, tmpdir):
        fname = tmpdir.join(self.net.project.name)
        f = op.io.project_to_hdf5(project=self.net.project, filename=fname)
        op.io.print_hdf5(f)
        filename = f.filename
        f.close()
        os.remove(filename)


if __name__ == '__main__':
    import py
    # All the tests in this file can be run with 'playing' this file
    t = HDF5Test()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
