import os
import py
import numpy as np
import openpnm as op
from openpnm.models.misc import from_neighbor_pores


class SalomeTest:

    def setup_class(self):
        np.random.seed(10)
        self.net = op.network.Cubic(shape=[2, 2, 2])
        self.net["pore.diameter"] = 0.5 + np.random.rand(self.net.Np) * 0.5
        Dt = from_neighbor_pores(self.net, prop="pore.diameter") * 0.5
        self.net["throat.diameter"] = Dt
        self.net["throat.length"] = 1.0

    def teardown_class(self):
        os.remove(f"{self.net.name}.py")
        os.remove("salome_custom.py")

    def test_export_data_salome(self):
        op.io.network_to_salome(network=self.net)
        assert os.path.isfile(f"{self.net.name}.py")
        op.io.network_to_salome(network=self.net, filename="salome_custom")
        assert os.path.isfile("salome_custom.py")


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = SalomeTest()
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
