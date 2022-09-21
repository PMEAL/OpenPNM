import os
import py
import pytest
import numpy as np
import openpnm as op
from openpnm.models.misc import from_neighbor_pores


@pytest.mark.skipif(os.name == 'nt', reason="Skipping on Windows")
class STLTest:

    def setup_class(self):
        np.random.seed(10)
        self.net = op.network.Cubic(shape=[2, 2, 2])
        self.net["pore.diameter"] = 0.5 + np.random.rand(self.net.Np) * 0.5
        Dt = from_neighbor_pores(self.net, prop="pore.diameter") * 0.5
        self.net["throat.diameter"] = Dt
        self.net["throat.length"] = 1.0

    def teardown_class(self):
        os.remove(f"{self.net.name}.stl")
        os.remove("custom_stl.stl")

    def test_export_data_stl(self):
        op.io.network_to_stl(network=self.net)
        assert os.path.isfile(f"{self.net.name}.stl")
        op.io.network_to_stl(network=self.net, filename="custom_stl")
        assert os.path.isfile("custom_stl.stl")


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = STLTest()
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
