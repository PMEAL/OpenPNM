import numpy as np
import openpnm as op
from openpnm.models import collections


class BCTest:
    def setup_class(self):
        np.random.seed(0)
        self.pn = op.network.Cubic(shape=[3, 3, 1], spacing=1e-4)
        self.pn.add_model_collection(collections.geometry.cones_and_cylinders())
        self.pn.regenerate_models()
        self.air = op.phase.Air(network=self.pn, name="air")
        self.air.add_model_collection(collections.physics.standard())
        self.air.regenerate_models()

    def test_add(self):
        fd = op.algorithms.FickianDiffusion(network=self.pn, phase=self.air)
        fd['pore.bc.rate'][1] = 1.0
        fd['pore.bc.value'][0] = 1.0
        fd.set_value_BC(pores=[0, 1, 2], values=3.0, mode='add')
        mask = np.isfinite(fd['pore.bc.value'])
        assert mask.sum() == 2
        assert fd['pore.bc.value'][mask].sum() == 6.0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 1
        fd.set_rate_BC(pores=[0, 1, 2], mode='remove')
        fd.set_value_BC(pores=[0, 1, 2], values=3.0, mode='add')
        mask = np.isfinite(fd['pore.bc.value'])
        assert mask.sum() == 3
        assert fd['pore.bc.value'][mask].sum() == 9.0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 0

    def test_overwrite(self):
        fd = op.algorithms.FickianDiffusion(network=self.pn, phase=self.air)
        fd['pore.bc.rate'][1] = 1.0
        fd['pore.bc.value'][0] = 1.0
        fd.set_value_BC(pores=[1, 2], values=3.0, mode='overwrite')
        mask = np.isfinite(fd['pore.bc.value'])
        assert mask.sum() == 3
        assert fd['pore.bc.value'][mask].sum() == 7.0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 0
        fd.set_rate_BC(pores=[0, 1, 2], rates=5.0, mode='overwrite')
        mask = np.isfinite(fd['pore.bc.rate'])
        assert mask.sum() == 3
        assert fd['pore.bc.rate'][mask].sum() == 15.0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 3
        assert np.isfinite(fd['pore.bc.value']).sum() == 0
        # Mimic removal of all bcs
        fd['pore.bc.rate'][1] = 1.0
        fd['pore.bc.value'][0] = 1.0
        fd.set_BC(pores=None, bctype=[], bcvalues=np.nan, mode='overwrite')
        assert np.isfinite(fd['pore.bc.rate']).sum() == 0
        assert np.isfinite(fd['pore.bc.value']).sum() == 0

    def test_remove(self):
        fd = op.algorithms.FickianDiffusion(network=self.pn, phase=self.air)
        # Remove given type from given pores
        fd['pore.bc.rate'][1] = 1.0
        fd['pore.bc.value'][0] = 1.0
        fd.set_value_BC(pores=[0, 1], mode='remove')
        mask = np.isfinite(fd['pore.bc.value'])
        assert mask.sum() == 0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 1
        # Remove given type from all pores
        fd['pore.bc.rate'][1] = 1.0
        fd['pore.bc.value'][0] = 1.0
        fd.set_BC(pores=None, bctype='rate', mode='remove')
        assert np.isfinite(fd['pore.bc.value']).sum() == 1
        assert np.isfinite(fd['pore.bc.rate']).sum() == 0
        # Remove all types from all pores
        fd['pore.bc.rate'][1] = 1.0
        fd['pore.bc.value'][0] = 1.0
        fd.set_BC(bctype=fd['pore.bc'].keys(), mode='remove')
        assert np.isfinite(fd['pore.bc.value']).sum() == 0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 0
        # Remove all types from all pores given bctype is []
        fd['pore.bc.rate'][1] = 1.0
        fd['pore.bc.value'][0] = 1.0
        fd.set_BC(mode='remove')
        assert np.isfinite(fd['pore.bc.value']).sum() == 0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 0
        # Remove all types from given pores
        fd['pore.bc.rate'][1] = 1.0
        fd['pore.bc.value'][0] = 1.0
        fd.set_BC(pores=[0, 1], mode='remove')
        assert np.isfinite(fd['pore.bc.value']).sum() == 0
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
        ad.set_outflow_BC(pores=self.pn.pores('back'), mode='overwrite')
        ad.run()

    def test_inlets_and_outlets(self):
        nwp = op.phase.GenericPhase(network=self.pn)
        nwp['throat.surface_tension'] = 0.480
        nwp['throat.contact_angle'] = 140
        nwp.add_model(propname='throat.entry_pressure',
                      model=op.models.physics.capillary_pressure.washburn)
        nwp.add_model(propname='pore.entry_pressure',
                      model=op.models.physics.capillary_pressure.washburn,
                      contact_angle=140,
                      surface_tension=0.480,
                      diameter='pore.diameter')

        drn = op.algorithms.Drainage(network=self.pn, phase=nwp)
        drn.set_inlets(pores=self.pn.pores('left'))
        pressures = np.logspace(np.log10(0.1e6), np.log10(8e6), 40)
        drn.run(pressures)
        drn.set_outlets(pores=self.pn.pores('right'))
        drn.apply_trapping()

    def test_modes_as_list(self):
        # check mode clear, with and without force
        fd = op.algorithms.FickianDiffusion(network=self.pn, phase=self.air)
        fd['pore.bc.rate'][1] = 1.0
        fd['pore.bc.value'][0] = 1.0
        fd.set_BC(bctype=['value', 'rate'], mode='clear')
        assert np.isfinite(fd['pore.bc.value']).sum() == 0
        assert np.isfinite(fd['pore.bc.rate']).sum() == 0


if __name__ == "__main__":

    t = BCTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith("test"):
            print("running test: " + item)
            t.__getattribute__(item)()
