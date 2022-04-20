import numpy as np
import scipy as sp
import openpnm as op
mgr = op.Workspace()


class PorosimetryTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[15, 15, 15], spacing=0.0005)
        self.geo = op.geometry.SpheresAndCylinders(network=self.net,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)
        self.hg = op.phase.Mercury(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.hg,
                                              geometry=self.geo)
        mod = op.models.physics.capillary_pressure.washburn
        self.phys.add_model(propname='throat.entry_pressure',
                            model=mod)

    # def test_no_late_filling(self):
    #     mip = op.metrics.Porosimetry(network=self.net, phase=self.hg)
    #     mip.set_inlets(pores=self.net.pores('left'))
    #     mip.run()
    #     assert len(np.unique(mip['pore.invasion_pressure'])) > 1
    #     assert len(np.unique(mip['pore.invasion_sequence'])) > 1
    #     assert len(np.unique(mip['throat.invasion_pressure'])) > 1
    #     assert len(np.unique(mip['throat.invasion_sequence'])) > 1

    # def test_late_pore_and_throat_filling(self):
    #     mip = op.metrics.Porosimetry(network=self.net, phase=self.hg)
    #     mip.set_inlets(pores=self.net.pores('left'))
    #     # Run without late pore filling
    #     mip.run()
    #     data_no_lpf = mip.get_intrusion_data()
    #     # Now run with late pore filling
    #     self.phys['pore.pc_star'] = 2/self.net['pore.diameter']
    #     self.phys.add_model(propname='pore.partial_filling',
    #                         pressure='pore.pressure',
    #                         Pc_star='pore.pc_star',
    #                         model=op.models.physics.multiphase.late_filling)
    #     mip.reset()
    #     mip.set_inlets(pores=self.net.pores('left'))
    #     mip.set_partial_filling(propname='pore.partial_filling')
    #     mip.run()
    #     self.phys.regenerate_models()
    #     data_w_lpf = mip.get_intrusion_data()
    #     assert np.all(np.array(data_w_lpf.Snwp) <= np.array(data_no_lpf.Snwp))
    #     # Now run with late throat filling
    #     self.phys['throat.pc_star'] = 2/self.net['throat.diameter']
    #     self.phys.add_model(propname='throat.partial_filling',
    #                         pressure='throat.pressure',
    #                         Pc_star='throat.pc_star',
    #                         model=op.models.physics.multiphase.late_filling)
    #     mip.reset()
    #     mip.set_inlets(pores=self.net.pores('left'))
    #     mip.set_partial_filling(propname='throat.partial_filling')
    #     mip.run()
    #     data_w_ltf = mip.get_intrusion_data()
    #     assert np.any(np.array(data_w_ltf.Snwp) <= np.array(data_w_lpf.Snwp))


if __name__ == '__main__':

    t = PorosimetryTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
