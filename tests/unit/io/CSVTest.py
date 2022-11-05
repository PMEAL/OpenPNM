import os
import openpnm as op


class CSVTest:

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
        self.phase_1['pore.baz'] = 11
        self.phase_1['throat.baz'] = 11
        self.phase_2['pore.baz'] = 12
        self.phase_2['throat.baz'] = 12

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_project_to_csv(self, tmpdir):
        fname = tmpdir.join(self.net.project.name)
        op.io.project_to_csv(self.net.project, filename=fname)
        assert os.path.isfile(fname + '.csv')
        os.remove(fname.dirpath().join(self.net.project.name + '.csv'))

    def test_network_to_csv(self, tmpdir):
        fname = tmpdir.join(self.net.name)
        op.io.network_to_csv(self.net, filename=fname)
        assert os.path.isfile(fname + '.csv')
        os.remove(fname.dirpath().join(self.net.name + '.csv'))

    def test_network_from_csv(self, tmpdir):
        fname = tmpdir.join(self.net.name)
        op.io.network_to_csv(self.net, filename=fname)
        net = op.io.network_from_csv(fname + '.csv')
        assert net.name == self.net.name
        assert net.Np == self.net.Np
        assert net.Nt == self.net.Nt
        assert net['pore.boo'].all() == self.net['pore.boo'].all()
        assert net['throat.boo'].all() == self.net['throat.boo'].all()
        os.remove(fname.dirpath().join(self.net.name + '.csv'))


if __name__ == '__main__':
    import py
    # All the tests in this file can be run with 'playing' this file
    t = CSVTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
