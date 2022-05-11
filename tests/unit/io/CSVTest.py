import openpnm as op
import pytest
import py
import os


class CSVTest:

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
        self.phase_1['pore.baz'] = 11
        self.phase_1['throat.baz'] = 11
        self.phase_2['pore.baz'] = 12
        self.phase_2['throat.baz'] = 12

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_save(self, tmpdir):
        fname = tmpdir.join(self.net.project.name)
        len_before = len(tmpdir.listdir())
        op.io.to_csv(network=self.net, phases=self.phase_1, filename=fname)
        assert len(tmpdir.listdir()) == (len_before + 1)
        os.remove(fname.dirpath().join(self.net.project.name + '.csv'))


if __name__ == '__main__':
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
