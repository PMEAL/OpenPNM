import openpnm as op
import scipy as sp
import matplotlib.pyplot as plt
mgr = op.Workspace()


class IPTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[10, 10, 10], spacing=0.0005)
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)
        self.water = op.phases.Water(network=self.net)
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.water,
                                              geometry=self.geo)
        mod = op.models.physics.capillary_pressure.washburn
        self.phys.add_model(propname='throat.entry_pressure',
                            model=mod)

    def test_set_inlets_overwrite(self):
        alg = op.algorithms.InvasionPercolation(network=self.net)
        alg.setup(phase=self.water)
        alg.set_inlets(pores=self.net.pores('top'))
        assert sp.sum(alg['pore.invasion_sequence'] == 0) == 100

        alg.set_inlets(pores=self.net.pores('bottom'))
        assert sp.sum(alg['pore.invasion_sequence'] == 0) == 200

        alg.set_inlets(pores=self.net.pores('top'), overwrite=True)
        assert sp.sum(alg['pore.invasion_sequence'] == 0) == 100

        alg.set_inlets(overwrite=True)
        assert sp.sum(alg['pore.invasion_sequence'] == 0) == 0

    def test_run(self):
        alg = op.algorithms.InvasionPercolation(network=self.net)
        alg.setup(phase=self.water)
        alg.set_inlets(pores=self.net.pores('top'))
        alg.run()
        assert alg['throat.invasion_sequence'].max() == (alg.Nt-1)

    def test_results(self):
        alg = op.algorithms.InvasionPercolation(network=self.net)
        alg.setup(phase=self.water)
        alg.set_inlets(pores=self.net.pores('top'))
        alg.run()
        d = alg.results(Snwp=0.5)
        assert set(d.keys()) == set(['pore.occupancy', 'throat.occupancy'])
        Vp = self.net['pore.volume']
        Vt = self.net['throat.volume']
        Vtot = Vp.sum() + Vt.sum()
        Vp_inv = Vp[d['pore.occupancy']].sum()
        Vt_inv = Vt[d['throat.occupancy']].sum()
        S = (Vp_inv + Vt_inv)/(Vtot)
        # Ensure saturation is close to requested value
        assert S < 0.6
        assert S > 0.4

    def test_trapping(self):
        alg = op.algorithms.InvasionPercolation(network=self.net)
        alg.setup(phase=self.water)
        alg.set_inlets(pores=self.net.pores('top'))
        alg.run()
        alg.apply_trapping(outlets=self.net.pores('bottom'))
        assert 'pore.trapped' in alg.labels()

    def test_plot_intrusion_curve(self):
        alg = op.algorithms.InvasionPercolation(network=self.net)
        alg.setup(phase=self.water)
        alg.set_inlets(pores=self.net.pores('top'))
        fig1 = alg.plot_intrusion_curve()
        assert fig1 is None
        alg.run()
        fig2 = alg.plot_intrusion_curve()
        ax2 = plt.gca()
        assert fig2 is not None
        alg.apply_trapping(outlets=self.net.pores('bottom'))
        fig3 = alg.plot_intrusion_curve()
        assert fig3 is not None
        ax3 = plt.gca()
        ydata2 = ax2.lines[0].get_ydata()
        ydata3 = ax3.lines[0].get_ydata()
        assert sp.any(ydata2-ydata3 != 0.0)
        plt.close('all')

if __name__ == '__main__':

    t = IPTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
