import pytest
import numpy as np
from numpy.testing import assert_approx_equal
import openpnm as op
import matplotlib.pyplot as plt


class IPTest:
    def setup_class(self):
        np.random.seed(0)
        self.net = op.network.Cubic(shape=[10, 10, 10], spacing=0.0005)
        self.net.add_model_collection(
            op.models.collections.geometry.spheres_and_cylinders)
        self.net.regenerate_models()
        self.water = op.phase.Water(network=self.net)
        self.water.add_model_collection(op.models.collections.physics.basic)
        self.water.regenerate_models()
        mod = op.models.physics.capillary_pressure.washburn
        self.water.add_model(propname="throat.entry_pressure", model=mod)

    def test_set_inlets(self):
        alg = op.algorithms.InvasionPercolation(network=self.net, phase=self.water)
        alg.set_inlet_BC(pores=self.net.pores("top"))
        assert np.sum(alg["pore.invasion_sequence"] == 0) == 100

        alg.set_inlet_BC(pores=self.net.pores("bottom"), mode='add')
        assert np.sum(alg["pore.invasion_sequence"] == 0) == 200

        alg.set_inlet_BC(mode='remove')
        assert np.sum(alg["pore.bc.inlet"] == True) == 0

        alg.set_inlet_BC(pores=self.net.pores("top"), mode='add')
        with pytest.raises(Exception):
            alg.set_outlet_BC(pores=self.net.pores("top"), mode='add')
        assert np.sum(alg["pore.bc.outlet"] == True) == 0
        alg.set_outlet_BC(pores=self.net.pores("top"), mode='overwrite')
        assert np.sum(alg["pore.bc.outlet"] == True) == 100

    def test_run(self):
        alg = op.algorithms.InvasionPercolation(network=self.net, phase=self.water)
        alg.set_inlet_BC(pores=self.net.pores("top"))
        alg.run()
        assert alg["throat.invasion_sequence"].max() == alg.Nt

    # def test_multiple_calls_to_run(self):
    #     alg = op.algorithms.InvasionPercolation(network=self.net, phase=self.water)
    #     alg.set_inlet_BC(pores=self.net.pores("top"))
    #     alg.run(n_steps=10)
    #     assert alg['pore.invasion_sequence'].max() == 10
    #     alg.run(n_steps=10)
    #     assert alg['pore.invasion_sequence'].max() == 20

    def test_results(self):
        alg = op.algorithms.InvasionPercolation(network=self.net, phase=self.water)
        alg.set_inlet_BC(pores=self.net.pores("top"))
        alg.run()
        assert alg['pore.invasion_sequence'].max() == 2695
        assert alg['pore.invasion_sequence'].min() == 0
        assert_approx_equal(alg['pore.invasion_pressure'].max(), 1967.22314587)
        assert alg['pore.invasion_pressure'].min() == 0

    def test_trapping(self):
        alg = op.algorithms.InvasionPercolation(network=self.net, phase=self.water)
        alg.set_inlet_BC(pores=self.net.pores("top"))
        alg.run()
        alg.set_outlet_BC(pores=self.net.pores("bottom"))
        alg.apply_trapping()
        assert "pore.trapped" in alg.keys()

    def test_plot_pc_curve(self):
        alg = op.algorithms.InvasionPercolation(network=self.net, phase=self.water)
        alg.set_inlet_BC(pores=self.net.pores("top"))
        with pytest.raises(KeyError):
            alg.pc_curve()
        alg.run()
        pc = alg.pc_curve()
        assert len(pc) == 2
        assert len(pc.pc) == 3700
        assert_approx_equal(pc.snwp.max(), 1.000000)


if __name__ == "__main__":

    t = IPTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith("test"):
            print("running test: " + item)
            t.__getattribute__(item)()
