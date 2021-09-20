import openpnm as op
import numpy as np


class BoundaryTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.net.add_boundary_pores()
        Ps_int = self.net.pores(labels=['top_boundary', 'bottom_boundary'],
                                mode='not')
        Ps_boun = self.net.pores(labels=['top_boundary', 'bottom_boundary'],
                                 mode='union')
        Pb_mask = np.random.random(len(Ps_boun)) < 0.5
        Ts_int = self.net.throats(labels=['internal', 'surface'], mode='union')
        TB_1 = self.net.find_neighbor_throats(pores=Ps_boun[Pb_mask])
        TB_2 = self.net.find_neighbor_throats(pores=Ps_boun[~Pb_mask])
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=Ps_int,
                                            throats=Ts_int)
        self.boun1 = op.geometry.Boundary(network=self.net,
                                          pores=Ps_boun[Pb_mask],
                                          throats=TB_1)
        self.boun2 = op.geometry.Boundary(network=self.net,
                                          pores=Ps_boun[~Pb_mask],
                                          throats=TB_2)
        self.geo.regenerate_models()
        self.boun1.regenerate_models()
        self.boun2.regenerate_models()

    def teardown_class(self):
        mgr = op.Workspace()
        mgr.clear()

    def test_plot_histogram(self):
        for obj in [self.boun1, self.boun2]:
            obj.show_hist()
            obj.show_hist(props=['pore.diameter', 'pore.volume',
                                 'throat.length'])
            obj.show_hist(props=['pore.diameter', 'pore.volume',
                                 'throat.length', 'throat.diameter',
                                 'pore.seed'])

    def test_boundary_with_alg(self):
        pn = op.network.Cubic(shape=[5, 5, 5], spacing=2.5e-5)
        pn.add_boundary_pores()
        Ps_int = pn.pores(labels=['*boundary'], mode='not')
        Ps_boun = pn.pores(labels=['*boundary'])
        Ts_int = pn.throats(labels=['*boundary'], mode='not')
        Ts_boun = pn.throats(labels=['*boundary'])
        geo = op.geometry.StickAndBall(network=pn,
                                       pores=Ps_int, throats=Ts_int)
        boun = op.geometry.Boundary(network=pn, pores=Ps_boun,
                                    throats=Ts_boun)
        geo.regenerate_models()
        boun.regenerate_models()
        air = op.phases.Air(network=pn)
        odiff = op.models.physics.diffusive_conductance.ordinary_diffusion
        phys_air_geo = op.physics.Standard(network=pn,
                                           phase=air,
                                           geometry=geo)
        phys_air_geo.add_model(propname="throat.diffusive_conductance",
                               model=odiff)
        phys_air_boun = op.physics.Standard(network=pn,
                                            phase=air,
                                            geometry=boun)
        phys_air_boun.add_model(propname="throat.diffusive_conductance",
                                model=odiff)
        phys_air_boun.regenerate_models()
        phys_air_geo.regenerate_models()
        health = phys_air_geo.check_data_health()
        for check in health.values():
            assert len(check) == 0
        checks = phys_air_boun.check_data_health().values()
        for check in checks:
            assert len(check) == 0
        FD = op.algorithms.FickianDiffusion(network=pn, phase=air)
        FD.set_value_BC(pores=pn.pores('top_boundary'), values=1.0)
        FD.set_value_BC(pores=pn.pores('bottom_boundary'), values=0.0)
        FD.run()
        assert FD.calc_effective_diffusivity() > 0
        SF = op.algorithms.StokesFlow(network=pn, phase=air)
        SF.set_value_BC(pores=pn.pores('top_boundary'), values=1.0)
        SF.set_value_BC(pores=pn.pores('bottom_boundary'), values=0.0)
        SF.run()
        assert SF.calc_effective_permeability() > 0


if __name__ == '__main__':

    t = BoundaryTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
