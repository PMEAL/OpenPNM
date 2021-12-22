import numpy as np
import numpy.testing as nt
import openpnm as op
import openpnm.models.geometry.diffusive_size_factors as gd
import openpnm.models.physics as pm


class TransientMultiPhysicsTest:

    def setup_class(self):
        np.random.seed(10)
        spacing=1e-1 
        # 2d network
        self.net = op.network.Cubic(shape=[10, 10, 1], spacing=spacing)
        # geometry
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo.add_model(propname='pore.diameter',
                           model=op.models.misc.constant,
                           value=spacing)
        self.geo.add_model(propname='throat.length',
                           model=op.models.misc.constant,
                           value=1e-15)
        self.geo.add_model(propname='throat.diameter',
                           model=op.models.misc.constant,
                           value=spacing)
        self.geo.add_model(propname='pore.cross_sectional_area',
                           model=op.models.misc.constant,
                           value=spacing)
        self.geo.add_model(propname='throat.cross_sectional_area',
                           model=op.models.misc.constant,
                           value=spacing)
        self.geo.add_model(propname='pore.volume',
                           model=op.models.misc.constant,
                           value=spacing**2)
        self.geo.add_model(propname='throat.diffusive_size_factors',
                           model=gd.squares_and_rectangles,
                           pore_diameter="pore.diameter",
                           throat_diameter="throat.diameter")
        # phase and physics
        self.air = op.phase.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net, 
                                              phase=self.air, 
                                              geometry=self.geo)
        # add diffusive conductance model
        self.air['pore.temperature'] = 300
        self.air.add_model(propname='pore.diffusivity',
                           model=op.models.misc.linear, 
                           m=1.860793056e-06,
                           b=-0.0005375624384,
                           prop='pore.temperature')
        self.phys.add_model(propname='throat.diffusive_conductance', 
                            model=pm.diffusive_conductance.generic_diffusive)
        # add thermal conductance model
        self.air.remove_model(propname='pore.thermal_conductivity')
        self.air.add_model(propname='pore.thermal_conductivity',
                           model=op.models.misc.constant,
                           value=0.0262,
                           regen_mode='constant')
        self.phys.add_model(propname='throat.thermal_conductance',
                            model=pm.thermal_conductance.generic_thermal) 
        # settings
        self.tfd_settings = {
            "conductance": "throat.diffusive_conductance",
            "quantity": "pore.concentration",
            "cache_A": False,
            "cache_b": False  
        }
        self.tfc_settings = {
            "conductance": "throat.thermal_conductance",
            "quantity": "pore.temperature",
            "pore_volume": "pore.heat_capacity",
            "cache_A": False,
            "cache_b": False  
        }
        self.pardiso = op.solvers.PardisoSpsolve()
        self.rk45 = op.integrators.ScipyRK45(verbose=True)
        
        # First algorithm, transient fourier conduction
        self.tfc = op.algorithms.TransientReactiveTransport(network=self.net,
                                                            phase=self.air)
        self.geo['pore.heat_capacity'] = self.geo['pore.volume']*1229.2875
        self.tfc.settings._update(self.tfc_settings)
        self.tfc.set_value_BC(self.net.pores("left"), 400)
        self.T0 = np.ones(self.tfc.Np) * 300
        self.T0[self.net.pores('left')] = 400

        # Second algorithm, transient fickian diffusion
        self.tfd = op.algorithms.TransientReactiveTransport(network=self.net,
                                                            phase=self.air)
        self.tfd.settings._update(self.tfd_settings)
        self.tfd.set_value_BC(self.net.pores("left"), 100)
        self.c0 = np.ones(self.tfd.Np) * 50
        self.c0[self.net.pores('left')] = 100

        # add variable props to algs
        self.tfd.set_variable_props('pore.temperature')
        self.tfc.set_variable_props('pore.concentration')
        
        # multiphysics algorithm
        algs = [self.tfc, self.tfd]
        self.tmp = op.algorithms.TransientMultiPhysics(algs, network=self.net)
        
    def test_run_algs(self):
        t_initial = 0
        t_final = 200
        t_step = 20
        n_steps = int((t_final - t_initial)/t_step) + 1
        t = np.linspace(t_initial, t_final, n_steps)
        y0 = np.hstack((self.T0, self.c0)) # ICs must include boundary condition
        tspan = [0, t_final]
        self.tmp.run(y0, tspan, saveat=t)
    
    def test_concentration(self):
        C = self.tmp.soln[self.net.Np:2*self.net.Np, :]
        C_avg = C[self.net.pores("left", mode='nor')].mean(axis=0)
        actual = C_avg
        desired = [50.00000000, 50.46652183, 51.00059257, 51.56467322, 
                   52.13159943, 52.68356036, 53.21031621, 53.70729529, 
                   54.17372048, 54.61117656, 55.02246821]
        nt.assert_allclose(actual, desired, rtol=1e-5)
        
if __name__ == '__main__':

    t = TransientMultiPhysicsTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
