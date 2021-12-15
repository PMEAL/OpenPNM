import openpnm as op
from pathlib import Path
import pickle
import os
import py


class PoreSpyTest:

    def setup_class(self):
        path = Path(os.path.realpath(__file__),
            '../../../../tests/fixtures/berea.net')
        with open(path, 'rb') as f:
            self.net = pickle.load(f)

    def test_load_PoreSpy_from_pickle(self, tmpdir):
        proj = op.io.from_porespy(self.net)
        net = proj.network
        assert net.Np == 1637
        assert net.Nt == 2785

    def test_load_PoreSpy_from_file(self, tmpdir):
        proj = op.io.from_porespy(filename=Path(os.path.realpath(__file__),
            '../../../../tests/fixtures/berea.net'))


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = PoreSpyTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
