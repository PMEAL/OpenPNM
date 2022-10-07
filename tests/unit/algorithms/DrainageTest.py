import pytest
import numpy as np
import openpnm as op
import matplotlib.pyplot as plt


class DrainageTest:
    def setup_class(self):
        self.pn = op.network.Demo(shape=[10, 10, 1], spacing=1e-4)
        self.air = op.phase.Air(network=self.pn)
        self.air.add_model_collection(op.models.collections.physics.standard)
        self.air['throat.surface_tension'] = 0.072
        self.air['throat.contact_angle'] = 160.0
        self.air.regenerate_models()

    def test_run(self):
        drn = op.algorithms.Drainage(network=self.pn, phase=self.air)
        drn.set_inlet_BC(pores=self.pn.pores('left'), mode='add')
        drn.run()
        assert np.all(drn['pore.invasion_pressure'][self.pn.conns].T
                      <= drn['throat.invasion_pressure'])
        # plt.imshow((drn['pore.invasion_pressure'] +
        #             20000*self.pn['pore.left']).reshape([10, 10]), origin='lower')

    def test_pccurve(self):
        drn = op.algorithms.Drainage(network=self.pn, phase=self.air)
        drn.set_inlet_BC(pores=self.pn.pores('left'), mode='add')
        drn.run()
        data = drn.pc_curve()
        assert np.all(np.unique(data[0])
                      == np.unique(drn['pore.invasion_pressure']))
        data = drn.pc_curve(np.linspace(0, 50000, 10))
        assert len(data[0]) == 10
        assert max(data[1]) == 1.0
        data = drn.pc_curve(np.linspace(0, 5000, 10))
        assert max(data[1]) < 1.0

    def test_apply_trapping(self):
        drn = op.algorithms.Drainage(network=self.pn, phase=self.air)
        drn.set_inlet_BC(pores=self.pn.pores('left'), mode='add')
        drn.set_outlet_BC(pores=self.pn.pores('right'), mode='add')
        drn.run()
        assert np.sum(np.isinf(drn['pore.invasion_pressure'])) > 0
        # plt.imshow((drn['pore.invasion_pressure'] +
        #             20000*self.pn['pore.left']).reshape([10, 10]), origin='lower')
        data = drn.pc_curve(np.linspace(0, 50000, 10))
        assert max(data[1]) < 1.0


if __name__ == "__main__":

    t = DrainageTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith("test"):
            print("running test: " + item)
            t.__getattribute__(item)()
