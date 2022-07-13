import openpnm as op
from numpy.testing import assert_allclose


class TransientAdvectionDiffusionTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[4, 3, 1], spacing=1.0)
        self.net.add_model_collection(
            op.models.collections.geometry.spheres_and_cylinders()
        )
        self.net.regenerate_models()
        self.phase = op.phase.Phase(network=self.net)
        self.phase['throat.diffusive_conductance'] = 1e-15
        self.phase['throat.hydraulic_conductance'] = 1e-15
        self.net['pore.volume'] = 1e-14
        self.sf = op.algorithms.StokesFlow(network=self.net, phase=self.phase)
        self.sf.settings._update({'quantity': 'pore.pressure',
                                  'conductance': 'throat.hydraulic_conductance'})
        self.sf.set_value_BC(pores=self.net.pores('right'), values=1)
        self.sf.set_value_BC(pores=self.net.pores('left'), values=0)
        self.sf.run()
        self.phase[self.sf.settings['quantity']] = self.sf.x
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phase.add_model(propname='throat.ad_dif_conductance', model=mod,
                             s_scheme='powerlaw')
        self.phase.regenerate_models()

    def test_transient_advection_diffusion(self):
        ad = op.algorithms.TransientAdvectionDiffusion(network=self.net,
                                                       phase=self.phase)
        ad.settings._update({
            'quantity': 'pore.concentration',
            'conductance': 'throat.ad_dif_conductance',
            'diffusive_conductance': 'throat.diffusive_conductance',
            'hydraulic_conductance': 'throat.hydraulic_conductance',
            'pressure': 'pore.pressure'
        })
        ad.set_value_BC(pores=self.net.pores('right'), values=2)
        ad.set_value_BC(pores=self.net.pores('left'), values=0)
        ad.run(x0=0, tspan=(0, 1))
        desired = 0.55642
        actual = ad.x.mean()
        assert_allclose(actual, desired, rtol=1e-5)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = TransientAdvectionDiffusionTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
