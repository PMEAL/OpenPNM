import openpnm as op
import scipy as sp


class HydraulicConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.diameter'] = 1.0
        self.geo['throat.diameter'] = 1.0
        self.geo['throat.length'] = 1.0e-9
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)

    def teardown_class(self):
        mgr = op.core.Workspace()
        mgr.clear()

    def test_hagen_poiseuille(self):
        mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
        self.phys.add_model(propname='throat.conductance1',
                            model=mod)
        assert sp.allclose(a=self.phys['throat.conductance1'][0],
                           b=1330.68207684)

        self.phys.add_model(propname='throat.conductance2',
                            model=mod,
                            calc_pore_len=True)
        assert sp.allclose(a=self.phys['throat.conductance2'][0],
                           b=1330.68207684)


if __name__ == '__main__':

    t = HydraulicConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
