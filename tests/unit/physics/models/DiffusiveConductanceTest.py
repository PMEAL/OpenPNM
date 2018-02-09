import openpnm
import openpnm.physics.models as gm
from numpy.testing import assert_approx_equal
from scipy.constants import pi


class DiffusiveConductanceTest:
    def setup_class(self):
        self.net = openpnm.network.Cubic(shape=[5, 5, 5], spacing=1e-8)
        self.geom = openpnm.geometry.GenericGeometry(network=self.net,
                                                     pores=self.net.Ps,
                                                     throats=self.net.Ts)
        self.geom['pore.diameter'] = 0.5e-8
        self.geom['pore.area'] = 0.25e-16
        self.geom['throat.diameter'] = 0.5e-8
        self.geom['throat.length'] = 0.5e-8
        self.geom['throat.area'] = pi * (0.5e-8)**2/4
        self.air = openpnm.phases.Air(network=self.net)
        self.phys = openpnm.physics.GenericPhysics(network=self.net,
                                                   phase=self.air,
                                                   geometry=self.geom)

    def test_bulk_diffusion(self):
        bulk_diffusion = gm.diffusive_conductance.bulk_diffusion
        self.phys.add_model(propname='throat.conductance',
                            model=bulk_diffusion)
        self.phys.regenerate_models(propnames='throat.conductance')
        assert_approx_equal(actual=self.phys['throat.conductance'][0],
                            desired=1.859713939028554e-12)

        self.phys.add_model(propname='throat.conductance',
                            model=bulk_diffusion,
                            calc_pore_len=True)
        self.phys.regenerate_models(propnames='throat.conductance')
        assert_approx_equal(actual=self.phys['throat.conductance'][0],
                            desired=1.859713939028554e-12)

    def test_mixed_diffusion(self):
        mixed_diffusion = gm.diffusive_conductance.mixed_diffusion
        self.phys.add_model(propname='throat.conductance',
                            model=mixed_diffusion)
        self.phys.regenerate_models(propnames='throat.conductance')
        assert_approx_equal(actual=self.phys['throat.conductance'][0],
                            desired=6.743639949710496e-14)
