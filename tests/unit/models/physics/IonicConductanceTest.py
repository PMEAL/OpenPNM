import openpnm as op
from openpnm.phases import mixtures
from numpy.testing import assert_allclose
import openpnm.models.geometry.conduit_lengths as _conduit_lengths

class IonicConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.diameter'] = 1.12
        self.geo['throat.diameter'] = 0.56
        self.geo['pore.area'] = 1
        self.geo['throat.area'] = 1
        self.geo['pore.volume']= 1
        self.geo['throat.volume']=1
        L1, Lt, L2 = _conduit_lengths.spheres_and_cylinders(self.geo).T
        self.geo['throat.conduit_lengths.pore1'] = L1
        self.geo['throat.conduit_lengths.throat'] = Lt
        self.geo['throat.conduit_lengths.pore2'] = L2
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.permittivity'] = 78
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        

    def test_ionic_conductance_poisson(self):
        # old poisson model, shape factor
        mpo = op.models.physics.poisson_shape_factors.ball_and_stick
        self.phys.add_model(propname="throat.poisson_shape_factors", model=mpo)
        mod1 = op.models.physics.ionic_conductance.poisson
        self.phys.add_model(propname='throat.ionic_conductance_from_mod',
                            model=mod1)
        self.phys.regenerate_models()
        # new poisson model, size factor
        mpo2 = op.models.geometry.diffusive_size_factors.spheres_and_cylinders
        self.geo.add_model(propname="throat.diffusive_size_factors",
                           model=mpo2)
        mod2 = op.models.physics.ionic_conductance.generic_ionic_poisson_laplace
        self.phys.add_model(propname='throat.ionic_conductance_generic',
                            model=mod2)
        self.phys.regenerate_models()
        actual = self.phys['throat.ionic_conductance_from_mod']
        desired = self.phys['throat.ionic_conductance_generic']
        assert_allclose(actual, desired, rtol=1e-5)
        
    def test_ionic_conductance_poisson_2D(self):
        self.geo['pore.area'] = self.geo['pore.diameter']
        self.geo['throat.area'] = self.geo['throat.diameter']
        # old poisson model, shape factor
        mpo = op.models.physics.poisson_shape_factors.ball_and_stick_2d
        self.phys.add_model(propname="throat.poisson_shape_factors", model=mpo)
        mod1 = op.models.physics.ionic_conductance.poisson
        self.phys.add_model(propname='throat.ionic_conductance_from_mod',
                            model=mod1)
        self.phys.regenerate_models()
        # new poisson model, size factor
        mpo2 = op.models.geometry.diffusive_size_factors.circles_and_rectangles
        self.geo.add_model(propname="throat.diffusive_size_factors",
                           model=mpo2)
        mod2 = op.models.physics.ionic_conductance.generic_ionic_poisson_laplace
        self.phys.add_model(propname='throat.ionic_conductance_generic',
                            model=mod2)
        self.phys.regenerate_models()
        actual = self.phys['throat.ionic_conductance_from_mod']
        desired = self.phys['throat.ionic_conductance_generic']
        assert_allclose(actual, desired, rtol=1e-5)

    def test_ionic_conductance_electroneutrality(self):
        # old electroneutrality model, shape factor
        sw = mixtures.SalineWater(network=self.net)
        Na = sw.components['Na_' + sw.name]
        Cl = sw.components['Cl_' + sw.name]
        H2O = sw.components['H2O_' + sw.name]
        # physics
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=sw, geometry=self.geo)
        sw['pore.concentration.'+Cl.name] = [10., 10., 10., 10., 10.,
                                             10., 10., 10., 10., 13.99930833,
                                             13.99930833, 13.99930833, 13.99930833, 13.99930833,
                                             13.99930833, 13.99930833, 13.99930833, 13.99930833,
                                             20., 20., 20., 20.,
                                             20., 20., 20., 20., 20.]

        sw['pore.concentration.'+Na.name] = [10., 10., 10., 10., 10.,
                                             10., 10., 10., 10., 13.9993082,
                                             13.9993082, 13.9993082, 13.9993082, 13.9993082, 13.9993082,
                                             13.9993082, 13.9993082, 13.9993082, 20., 20.,
                                             20., 20., 20., 20., 20.,
                                             20., 20.]
        mpo = op.models.physics.poisson_shape_factors.ball_and_stick
        self.phys.add_model(propname="throat.poisson_shape_factors", model=mpo)
        current = op.models.physics.ionic_conductance.electroneutrality
        self.phys.add_model(propname='throat.ionic_conductance',
                            ions=[Na.name, Cl.name],
                            model=current, regen_mode='normal')
        mpo2 = op.models.geometry.diffusive_size_factors.spheres_and_cylinders
        self.geo.add_model(propname="throat.diffusive_size_factors",
                           model=mpo2)
        current2 = op.models.physics.ionic_conductance.generic_ionic_electroneutrality
        self.phys.add_model(propname='throat.ionic_conductance_generic',
                            ions=[Na.name, Cl.name],
                            model=current2, regen_mode='normal')
        self.phys.regenerate_models()
        actual = self.phys['throat.ionic_conductance']
        desired = self.phys['throat.ionic_conductance_generic']
        assert_allclose(actual, desired, rtol=1e-5)


if __name__ == '__main__':

    t = IonicConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()