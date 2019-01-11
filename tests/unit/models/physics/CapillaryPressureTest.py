import openpnm as op
import scipy as sp
from numpy.testing import assert_approx_equal


class CapillaryPressureTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['throat.diameter'] = 1
        self.geo['pore.diameter'] = 1
        self.water = op.phases.GenericPhase(network=self.net)
        self.water['pore.surface_tension'] = 0.072
        self.water['pore.contact_angle'] = 120
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              geometry=self.geo,
                                              phase=self.water)

    def test_washburn_pore_values(self):
        f = op.models.physics.capillary_pressure.washburn
        self.phys.add_model(propname='pore.capillary_pressure',
                            model=f,
                            surface_tension='pore.surface_tension',
                            contact_angle='pore.contact_angle',
                            diameter='pore.diameter')
        self.phys.add_model(propname='throat.capillary_pressure',
                            model=f,
                            surface_tension='pore.surface_tension',
                            contact_angle='pore.contact_angle',
                            diameter='throat.diameter')
        self.phys.regenerate_models()
        a = 0.14399999999999993
        assert_approx_equal(self.water['pore.capillary_pressure'][0], a)
        self.phys.remove_model('pore.capillary_pressure')

    def test_washburn_throat_values(self):
        f = op.models.physics.capillary_pressure.washburn
        self.phys.add_model(propname='throat.capillary_pressure',
                            model=f,
                            surface_tension='throat.surface_tension',
                            contact_angle='throat.contact_angle',
                            diameter='throat.diameter')
        self.phys.regenerate_models()
        a = 0.14399999999999993
        assert_approx_equal(self.water['throat.capillary_pressure'][0], a)
        self.phys.remove_model('throat.capillary_pressure')

    def test_purcell_pore_values(self):
        f = op.models.physics.capillary_pressure.purcell
        self.phys.add_model(propname='pore.capillary_pressure',
                            model=f,
                            r_toroid=0.1,
                            surface_tension='pore.surface_tension',
                            contact_angle='pore.contact_angle',
                            diameter='pore.diameter')
        self.phys.regenerate_models()
        a = 0.2648436086476371
        assert_approx_equal(self.water['pore.capillary_pressure'][0], a)
        self.phys.remove_model('pore.capillary_pressure')

    def test_purcell_throat_values(self):
        f = op.models.physics.capillary_pressure.purcell
        self.phys.add_model(propname='throat.capillary_pressure',
                            model=f,
                            r_toroid=0.1,
                            surface_tension='pore.surface_tension',
                            contact_angle='pore.contact_angle',
                            diameter='throat.diameter')
        self.phys.regenerate_models()
        a = 0.2648436086476371
        assert_approx_equal(self.water['throat.capillary_pressure'][0], a)
        self.phys.remove_model('throat.capillary_pressure')

    def test_purcell_bidirectional(self):
        f = op.models.physics.capillary_pressure.purcell_bidirectional
        self.geo['pore.touch'] = (sp.random.random(self.geo.Np)+0.5)*0.1
        self.phys.add_model(propname='throat.bidirectional',
                            model=f,
                            r_toroid=0.1,
                            surface_tension='pore.surface_tension',
                            contact_angle='pore.contact_angle',
                            pore_diameter='pore.touch')
        diff = (self.phys['throat.bidirectional'][:, 0] -
                self.phys['throat.bidirectional'][:, 1])
        assert sp.any(diff != 0)

    def test_sinusoidal_bidirectional(self):
        f = op.models.physics.capillary_pressure.sinusoidal_bidirectional
        self.geo['pore.touch'] = (sp.random.random(self.geo.Np)+0.5)*0.1
        self.geo['throat.length'] = 1.0
        self.phys.add_model(propname='throat.bidirectional',
                            model=f,
                            r_toroid=0.25,
                            surface_tension='pore.surface_tension',
                            contact_angle='pore.contact_angle',
                            pore_diameter='pore.touch')
        diff = (self.phys['throat.bidirectional'][:, 0] -
                self.phys['throat.bidirectional'][:, 1])
        assert sp.any(diff != 0)

    def test_ransohoff_snapoff_verts(self):
        ws = op.Workspace()
        ws.clear()
        bp = sp.array([[0.25, 0.25, 0.25], [0.25, 0.75, 0.25],
                       [0.75, 0.25, 0.25], [0.75, 0.75, 0.25],
                       [0.25, 0.25, 0.75], [0.25, 0.75, 0.75],
                       [0.75, 0.25, 0.75], [0.75, 0.75, 0.75]])
        scale = 1e-4
        sp.random.seed(1)
        p = (sp.random.random([len(bp), 3])-0.5)/1000
        bp += p
        fiber_rad = 2e-6
        bp = op.topotools.reflect_base_points(bp, domain_size=[1, 1, 1])
        prj = op.materials.VoronoiFibers(fiber_rad=fiber_rad,
                                         resolution=1e-6,
                                         shape=[scale, scale, scale],
                                         points=bp*scale,
                                         name='test')
        net = prj.network
        del_geom = prj.geometries()['test_del']
        vor_geom = prj.geometries()['test_vor']
        f = op.models.physics.capillary_pressure.ransohoff_snap_off
        water = op.phases.GenericPhase(network=net)
        water['pore.surface_tension'] = 0.072
        water['pore.contact_angle'] = 45
        phys1 = op.physics.GenericPhysics(network=net,
                                          geometry=del_geom,
                                          phase=water)
        phys1.add_model(propname='throat.snap_off',
                        model=f,
                        wavelength=fiber_rad)
        phys1.add_model(propname='throat.snap_off_pair',
                        model=f,
                        wavelength=fiber_rad,
                        require_pair=True)
        phys2 = op.physics.GenericPhysics(network=net,
                                          geometry=vor_geom,
                                          phase=water)
        phys2.add_model(propname='throat.snap_off',
                        model=f,
                        wavelength=fiber_rad)
        phys2.add_model(propname='throat.snap_off_pair',
                        model=f,
                        wavelength=fiber_rad,
                        require_pair=True)
        ts = ~net['throat.interconnect']
        assert ~sp.any(sp.isnan(water['throat.snap_off'][ts]))
        assert sp.any(sp.isnan(water['throat.snap_off_pair'][ts]))
        assert sp.any(~sp.isnan(water['throat.snap_off_pair'][ts]))


if __name__ == '__main__':

    t = CapillaryPressureTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
