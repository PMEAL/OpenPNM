import openpnm
from numpy.testing import assert_approx_equal


class DiffusionTest:
    def setup_class(self):
        self.net = openpnm.network.Cubic(shape=[5, 5, 5], spacing=1e-8)
        self.geom = openpnm.geometry.GenericGeometry(network=self.net,
                                                     pores=self.net.Ps,
                                                     throats=self.net.Ts)
        self.geom['pore.diameter'] = 1e-8
        self.geom['pore.area'] = 1e-16
        self.geom['throat.diameter'] = 1e-8
        self.geom['throat.length'] = 1e-9
        self.geom['throat.area'] = 1e-16
        self.air = openpnm.phases.Air(network=self.net)
        self.phys = openpnm.physics.GenericPhysics(network=self.net,
                                                   phase=self.air,
                                                   geometry=self.geom)

    def test_knudsen(self):
        knudsen = openpnm.physics.models.diffusion.knudsen
        self.phys.add_model(propname='pore.knudsen_diffusivity',
                            model=knudsen)
        self.phys.regenerate_models(propnames='pore.knudsen_diffusivity')
        assert_approx_equal(actual=self.phys['pore.knudsen_diffusivity'][0],
                            desired=1.5558748999706673e-06)
#        print(self.air['pore.diffusivity'])

    def test_knudsen_scaling(self):
        knudsen_scaling = openpnm.physics.models.diffusion.knudsen_scaling
        self.phys.add_model(propname='pore.mixed_diffusivity',
                            model=knudsen_scaling)
        self.phys.regenerate_models(propnames='pore.mixed_diffusivity')
        print(self.phys['pore.mixed_diffusivity'][0])

        assert_approx_equal(actual=self.phys['pore.mixed_diffusivity'][0],
                            desired=1.4469860405433859e-06)
  