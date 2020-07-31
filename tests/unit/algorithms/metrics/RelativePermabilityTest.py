import openpnm as op
import numpy.testing as nt
mgr = op.Workspace()


class RelativePermeabilityTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1)
        self.geo = op.geometry.GenericGeometry(network=self.net)
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
        self.inlet_pores = self.net.pores('left')
        ip = op.algorithms.InvasionPercolation(network=self.net,
                                               phase=self.non_wet_phase)
        ip.set_inlets(pores=self.inlet_pores)
        ip.run()
        self.non_wet_phase.update(ip.results())

    def test_one_phase_definition(self):
        rp = op.algorithms.metrics.RelativePermeability(network=self.net)
        rp.setup(invading_phase=self.non_wet_phase,
                 invasion_sequence='invasion_sequence')
        rp.run(Snwp_num=10)
        results = rp.get_Kr_data()
        assert results['relperm_wp'] == {}

    def test_overwriting_boundary_faces(self):
        inlets = {'x': 'left', 'y': 'left', 'z': 'left'}
        outlets = {'x': 'right', 'y': 'right', 'z': 'right'}
        rp = op.algorithms.metrics.RelativePermeability(network=self.net)
        rp.setup(invading_phase=self.non_wet_phase,
                 defending_phase=self.wet_phase,
                 invasion_sequence='invasion_sequence',
                 flow_inlets=inlets,
                 flow_outlets=outlets)
        rp.run(Snwp_num=10)
        results = rp.get_Kr_data()
        kx = results['relperm_wp']['x']
        ky = results['relperm_wp']['y']
        kr = [0.7230822778535343, 0.5469031280514686, 0.4675498520331332,
              0.10041453914739418, 1.2428494917580884e-06, 1.0000000000000004e-06,
              1.0000000000000004e-06, 1.0000000000000004e-06, 1.0000000000000004e-06,
              1.0000000000000004e-06]
        nt.assert_allclose(kx, ky)
        nt.assert_allclose(kx, kr)

    def test_lacking_boundary_faces(self):
        inlets = {'x': 'top'}
        outlets = {'x': 'bottom'}
        rp = op.algorithms.metrics.RelativePermeability(network=self.net)
        rp.setup(invading_phase=self.non_wet_phase,
                 defending_phase=self.wet_phase,
                 invasion_sequence='invasion_sequence',
                 flow_inlets=inlets,
                 flow_outlets=outlets)
        rp.run(Snwp_num=10)
        results = rp.get_Kr_data()
        kx = results['kr_wp']['x']
        kz = results['kr_wp']['z']
        kr = [0.5953556221922877, 0.42713264157774794, 0.3658925423425995,
              0.21493111700350034, 1.2600781827032384e-06, 1.000000000000001e-06,
              1.000000000000001e-06, 1.000000000000001e-06, 1.000000000000001e-06,
              1.000000000000001e-06]
        nt.assert_allclose(kx, kz)
        nt.assert_allclose(kx, kr)

    def test_user_defined_boundary_face(self):
        pores_in = self.net.pores('top')
        pores_out = self.net.pores('bottom')
        self.net.set_label(pores=pores_in, label='pore_in')
        self.net.set_label(pores=pores_out, label='pore_out')
        inlets = {'x': 'pore_in'}
        outlets = {'x': 'pore_out'}
        rp = op.algorithms.metrics.RelativePermeability(network=self.net)
        rp.setup(invading_phase=self.non_wet_phase,
                  defending_phase=self.wet_phase,
                  invasion_sequence='invasion_sequence',
                  flow_inlets=inlets,
                  flow_outlets=outlets)
        rp.run(Snwp_num=10)
        results = rp.get_Kr_data()
        kx = results['kr_wp']['x']
        kz = results['kr_wp']['z']
        kr = [0.5953556221922877, 0.42713264157774794, 0.3658925423425995,
              0.21493111700350034, 1.2600781827032384e-06, 1.000000000000001e-06,
              1.000000000000001e-06, 1.000000000000001e-06, 1.000000000000001e-06,
              1.000000000000001e-06]
        nt.assert_allclose(kx, kz)
        nt.assert_allclose(kx, kr)


if __name__ == '__main__':

    t = RelativePermeabilityTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
