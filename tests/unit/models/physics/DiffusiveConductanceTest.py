import openpnm as op
import scipy as sp


class DiffusiveConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.diameter'] = 1.0
        self.geo['pore.area'] = 1.0
        self.geo['throat.diameter'] = 1.0
        self.geo['throat.length'] = 1e-9
        self.geo['throat.area'] = 1
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)

    def test_bulk_diffusion(self):
        mod = op.models.physics.diffusive_conductance.bulk_diffusion
        self.phys.add_model(propname='throat.conductance1',
                            model=mod)
        self.phys.regenerate_models()
        assert sp.allclose(a=self.phys['throat.conductance1'][0],
                           b=0.00084552)

        self.phys.add_model(propname='throat.conductance2',
                            model=mod,
                            calc_pore_len=True)
        self.phys.regenerate_models()
        assert sp.allclose(a=self.phys['throat.conductance2'][0],
                           b=0.00084552)


if __name__ == '__main__':

    t = DiffusiveConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
