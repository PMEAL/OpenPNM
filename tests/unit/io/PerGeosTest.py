import os
from pathlib import Path

import openpnm as op


class PerGeosTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.clear()

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_load_PerGeos_simple(self, tmpdir):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/PerGeos/simplePNM.am')
        network = op.io.network_from_pergeos(path)
        assert network.Np == 3
        assert network.Nt == 3

    def test_load_PerGeos_mandatory(self, tmpdir):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/PerGeos/mandatory.am')
        project = op.io.network_from_pergeos(path)
        network = project.network
        assert network.Np == 3
        assert network.Nt == 3

    def test_load_PerGeos_flooded(self, tmpdir):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/PerGeos/flooded.am')
        network = op.io.network_from_pergeos(path)
        assert network.Np == 225
        assert network.Nt == 301

    def test_save_PerGeos(self, tmpdir):
        net = op.network.Cubic(shape=[5, 5, 5])
        fname = tmpdir.join(net.project.name)
        len_before = len(tmpdir.listdir())
        op.io.network_to_pergeos(net, filename=fname)
        assert len(tmpdir.listdir()) == (len_before + 1)
        os.remove(fname.dirpath().join(net.project.name + '.am'))


if __name__ == '__main__':
    import py

    # All the tests in this file can be run with 'playing' this file
    t = PerGeosTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
