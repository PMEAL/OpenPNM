import numpy as np
import openpnm as op
from openpnm.models import collections


class BCTest:
    def setup_class(self):
        self.pn = op.network.Cubic(shape=[3, 3, 1], spacing=1e-4)
        self.pn.add_model_collection(collections.geometry.cones_and_cylinders())
        self.pn.regenerate_models()
        self.air = op.phase.Air(network=self.pn, name="air")
        self.air.add_model_collection(collections.physics.standard())
        self.air.regenerate_models()

    def test_add(self):
        # check mode add, with and without force
        fd = op.algorithms.FickianDiffusion(network=self.pn, phase=self.air)
        fd['pore.bc.rate'][1] = 1.0
        fd['pore.bc.value'][0] = 1.0
        fd.set_value_BC(pores=[0, 1, 2], values=3.0, mode='add')
        mask = np.isfinite(fd['pore.bc.value'])
        assert mask.sum() == 2
        assert fd['pore.bc.value'][mask].sum() == 4.0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 1
        fd.set_value_BC(pores=[0, 1, 2], values=3.0, mode='add', force=True)
        mask = np.isfinite(fd['pore.bc.value'])
        assert mask.sum() == 3
        assert fd['pore.bc.value'][mask].sum() == 7.0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 0

    def test_overwite(self):
        # check mode overwrite, with and without force
        fd = op.algorithms.FickianDiffusion(network=self.pn, phase=self.air)
        fd['pore.bc.rate'][1] = 1.0
        fd['pore.bc.value'][0] = 1.0
        fd.set_value_BC(pores=[0, 1, 2], values=3.0, mode='overwrite')
        mask = np.isfinite(fd['pore.bc.value'])
        assert mask.sum() == 2
        assert fd['pore.bc.value'][mask].sum() == 6.0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 1
        fd.set_value_BC(pores=[0, 1, 2], values=3.0, mode='overwrite', force=True)
        mask = np.isfinite(fd['pore.bc.value'])
        assert mask.sum() == 3
        assert fd['pore.bc.value'][mask].sum() == 9.0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 0

    def test_remove(self):
        # check mode remove, with and without force
        fd = op.algorithms.FickianDiffusion(network=self.pn, phase=self.air)
        fd['pore.bc.rate'][1] = 1.0
        fd['pore.bc.value'][0] = 1.0
        fd.set_value_BC(pores=[0, 1], mode='remove')
        mask = np.isfinite(fd['pore.bc.value'])
        assert mask.sum() == 0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 1
        fd.set_value_BC(pores=[0, 1], mode='remove', force=True)
        assert np.isfinite(fd['pore.bc.value']).sum() == 0
        assert np.isfinite(fd['pore.bc.value']).sum() == 0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 0

    def test_clear(self):
        # check mode clear, with and without force
        fd = op.algorithms.FickianDiffusion(network=self.pn, phase=self.air)
        fd['pore.bc.rate'][1] = 1.0
        fd['pore.bc.value'][0] = 1.0
        fd.set_value_BC(mode='clear')
        assert np.isfinite(fd['pore.bc.value']).sum() == 0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 1
        fd.set_value_BC(mode='clear', force=True)
        assert np.isfinite(fd['pore.bc.rate']).sum() == 0

    def test_outflow(self):
        del self.air.models['throat.diffusive_conductance']
        del self.air.models['throat.hydraulic_conductance']
        self.air['throat.diffusive_conductance'] = 1e-15
        self.air['throat.hydraulic_conductance'] = 1e-15
        flow = op.algorithms.StokesFlow(network=self.pn, phase=self.air)
        flow.set_value_BC(pores=self.pn.pores('left'), values=1)
        flow.set_value_BC(pores=self.pn.pores('right'), values=0)
        flow.run()
        ad = op.algorithms.AdvectionDiffusion(network=self.pn, phase=self.air)
        ad.settings['cache'] = False
        ad.set_value_BC(pores=self.pn.pores('front'), values=1)
        ad.set_outflow_BC(pores=self.pn.pores('back'), mode='overwrite', force=False)
        ad.run()


if __name__ == "__main__":

    t = BCTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith("test"):
            print("running test: " + item)
            t.__getattribute__(item)()
