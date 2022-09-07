import os
import py
import numpy as np
import openpnm as op
import pytest


class COMSOLTest:

    def setup_class(self):
        np.random.seed(10)
        self.net2d = op.network.Cubic(shape=[3, 4])
        self.net2d["pore.diameter"] = np.random.rand(self.net2d.Np)
        self.net2d["throat.diameter"] = np.random.rand(self.net2d.Nt)
        self.net3d = op.network.Cubic(shape=[3, 4, 5])
        self.net3d["pore.diameter"] = np.random.rand(self.net3d.Np)
        self.net3d["throat.diameter"] = np.random.rand(self.net3d.Nt)


    def teardown_class(self):
        os.remove(f"{self.net2d.name}.mphtxt")

    def test_export_data_2d_network(self):
        op.io.network_to_comsol(network=self.net2d)
        assert os.path.isfile(f"{self.net2d.name}.mphtxt")

    def test_export_data_3d_network(self):
        with pytest.raises(Exception):
            op.io.network_to_comsol(network=self.net3d)


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = COMSOLTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
    t.teardown_class()
