import openpnm as op
import pytest
import py
import os
from pathlib import Path


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
        project = op.io.from_pergeos(path)
        network = project.network
        assert network.Np == 3
        assert network.Nt == 3

    def test_load_PerGeos_mandatory(self, tmpdir):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/PerGeos/mandatory.am')
        project = op.io.from_pergeos(path)
        network = project.network
        assert network.Np == 3
        assert network.Nt == 3

    def test_load_PerGeos_flooded(self, tmpdir):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/PerGeos/flooded.am')
        project = op.io.from_pergeos(path)
        network = project.network
        assert network.Np == 225
        assert network.Nt == 301

    def test_save_PerGeos(self, tmpdir):
        net = op.network.Cubic(shape=[5, 5, 5])
        water = op.phase.Water(network=net)
        fname = tmpdir.join(net.project.name)
        len_before = len(tmpdir.listdir())
        op.io.to_pergeos(network=net, phases=water, filename=fname)
        print(tmpdir.listdir())
        print(len(tmpdir.listdir()))
        print(len_before)
        assert len(tmpdir.listdir()) == (len_before + 1)
        os.remove(fname.dirpath().join(net.project.name + '.am'))


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = PerGeosTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
