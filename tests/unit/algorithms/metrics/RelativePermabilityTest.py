import openpnm as op
import numpy.testing as nt
mgr = op.Workspace()


class RelativePermeabilityTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo["pore.diameter"] = 0.5
        self.geo["pore.area"] = 0.5**2
        self.geo["pore.volume"] = 0.5**3
        self.geo["throat.diameter"] = 0.3
        self.geo["throat.area"] = 0.3**2
        self.geo["throat.volume"] = 0.3**3
        self.geo["throat.conduit_lengths.throat"] = 1
        self.geo["throat.conduit_lengths.pore1"] = 0.3
        self.geo["throat.conduit_lengths.pore2"] = 0.9
        self.non_wet_phase = op.phases.Air(network=self.net)
        self.wet_phase = op.phases.Water(network=self.net)
        mod = op.models.geometry.hydraulic_size_factors.spheres_and_cylinders
        self.geo.add_model(propname='throat.hydraulic_size_factors', model=mod)
        mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
        self.non_wet_phase.add_model(propname='throat.hydraulic_conductance',
                                     model=mod)
        self.wet_phase.add_model(propname='throat.hydraulic_conductance',
                                 model=mod)
        mod = op.models.physics.capillary_pressure.washburn
        self.non_wet_phase.add_model(propname='throat.entry_pressure',
                                     model=mod)
        self.wet_phase.add_model(propname='throat.entry_pressure',
                                 model=mod)
        self.inlet_pores = self.net.pores('back')
        ip = op.algorithms.InvasionPercolation(network=self.net,
                                               phase=self.non_wet_phase)
        ip.set_inlets(pores=self.inlet_pores)
        ip.run()
        self.non_wet_phase.update(ip.results())

    def test_one_phase_definition(self):
        rp = op.algorithms.metrics.RelativePermeability(network=self.net)
        rp.settings.update({'nwp': self.non_wet_phase.name,
                            'invasion_sequence': 'invasion_sequence'})
        rp.run(Snwp_num=10)
        results = rp.get_Kr_data()
        assert results['kr_wp'] is None

    def test_overwriting_boundary_faces(self):
        inlets = {'x': 'back', 'y': 'back', 'z': 'back'}
        outlets = {'x': 'front', 'y': 'front', 'z': 'front'}
        rp = op.algorithms.metrics.RelativePermeability(network=self.net)
        rp.settings.update({'nwp': self.non_wet_phase.name,
                            'wp': self.wet_phase.name,
                            'invasion_sequence': 'invasion_sequence'
                            })
        rp.settings['flow_inlets'].update(inlets)
        rp.settings['flow_outlets'].update(outlets)
        rp.run(Snwp_num=10)
        results = rp.get_Kr_data()
        kx = results['kr_wp']['x']
        ky = results['kr_wp']['y']
        kr = [7.003833e-01, 4.675499e-01, 4.675499e-01, 2.371033e-06,
              1.216706e-06, 1.000000e-06, 1.000000e-06, 1.000000e-06,
              1.000000e-06, 1.000000e-06]
        nt.assert_allclose(kx, ky, rtol=1e-6)
        nt.assert_allclose(kx, kr, rtol=1e-6)

    def test_lacking_boundary_faces(self):
        rp = op.algorithms.metrics.RelativePermeability(network=self.net)
        inlets = {'x': 'top'}
        outlets = {'x': 'bottom'}
        rp.settings.update({'nwp': self.non_wet_phase.name,
                            'wp': self.wet_phase.name,
                            'invasion_sequence': 'invasion_sequence',
                            })
        rp.settings['flow_inlets'].update(inlets)
        rp.settings['flow_outlets'].update(outlets)
        rp.run(Snwp_num=10)
        results = rp.get_Kr_data()
        kx = results['kr_wp']['x']
        kz = results['kr_wp']['z']
        kr = [5.982845e-01, 4.060000e-01, 4.060000e-01, 2.046288e-01,
              1.065283e-06, 1.000000e-06, 1.000000e-06, 1.000000e-06,
              1.000000e-06, 1.000000e-06]
        nt.assert_allclose(kx, kz, rtol=1e-6)
        nt.assert_allclose(kx, kr, rtol=1e-6)

    def test_user_defined_boundary_face(self):
        pores_in = self.net.pores('top')
        pores_out = self.net.pores('bottom')
        self.net.set_label(pores=pores_in, label='pore_in')
        self.net.set_label(pores=pores_out, label='pore_out')
        inlets = {'x': 'pore_in'}
        outlets = {'x': 'pore_out'}
        rp = op.algorithms.metrics.RelativePermeability(network=self.net)
        rp.settings.update({'nwp': self.non_wet_phase.name,
                            'wp': self.wet_phase.name,
                            'invasion_sequence': 'invasion_sequence'
                            })
        rp.settings['flow_inlets'].update(inlets)
        rp.settings['flow_outlets'].update(outlets)
        rp.run(Snwp_num=10)
        results = rp.get_Kr_data()
        kx = results['kr_wp']['x']
        kz = results['kr_wp']['z']
        kr = [5.982845e-01, 4.060000e-01, 4.060000e-01, 2.046288e-01,
              1.065283e-06, 1.000000e-06, 1.000000e-06, 1.000000e-06,
              1.000000e-06, 1.000000e-06]
        nt.assert_allclose(kx, kz, rtol=1e-6)
        nt.assert_allclose(kx, kr, rtol=1e-6)

    def setup_model2d(self, shape):
        self.net = op.network.Cubic(shape=shape, spacing=0.0005)
        self.geo = op.geometry.SpheresAndCylinders(network=self.net,
                                                   pores=self.net.Ps,
                                                   throats=self.net.Ts)
        self.non_wet_phase = op.phases.Air(network=self.net)
        self.wet_phase = op.phases.Water(network=self.net)
        mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
        self.non_wet_phase.add_model(propname='throat.hydraulic_conductance',
                                     model=mod)
        self.wet_phase.add_model(propname='throat.hydraulic_conductance',
                                 model=mod)
        mod = op.models.physics.capillary_pressure.washburn
        self.non_wet_phase.add_model(propname='throat.entry_pressure',
                                     model=mod)
        self.wet_phase.add_model(propname='throat.entry_pressure',
                                 model=mod)
        if shape[1] != 1:
            self.inlet_pores = self.net.pores('back')
        else:
            self.inlet_pores = self.net.pores('left')
        ip = op.algorithms.InvasionPercolation(network=self.net,
                                               phase=self.non_wet_phase)
        ip.set_inlets(pores=self.inlet_pores)
        ip.run()
        self.non_wet_phase.update(ip.results())

    def test_model2d_one_phase_curve(self):
        for i in range(3):
            shape = [10, 10, 10]
            shape[i] = 1
            self.setup_model2d(shape=shape)
            rp = op.algorithms.metrics.RelativePermeability(network=self.net)
            rp.settings.update({'nwp': self.non_wet_phase.name,
                                'invasion_sequence': 'invasion_sequence'})
            rp.run(Snwp_num=10)
            results = rp.get_Kr_data()
            assert results['kr_wp'] is None
            nt.assert_allclose(len(results['sat']), 2)

    def test_model2d_two_phase_curve(self):
        for i in range(3):
            shape = [10, 10, 10]
            shape[i] = 1
            self.setup_model2d(shape=shape)
            rp = op.algorithms.metrics.RelativePermeability(network=self.net)
            rp.settings.update({'nwp': self.non_wet_phase.name,
                                'wp': self.wet_phase.name,
                                'invasion_sequence': 'invasion_sequence'})
            rp.run(Snwp_num=10)
            results = rp.get_Kr_data()
            nt.assert_allclose(len(results['kr_wp']), 2)
            nt.assert_allclose(len(results['kr_nwp']), 2)
            nt.assert_allclose(len(results['sat']), 2)


if __name__ == '__main__':

    t = RelativePermeabilityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
