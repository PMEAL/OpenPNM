import openpnm as op
import numpy as np
from openpnm.algorithms import MixedInvasionPercolationCoop as mpc
import matplotlib.pyplot as plt
import openpnm.models.geometry as gm


plt.close('all')
wrk = op.Workspace()
wrk.loglevel = 50


class MixedPercolationCoopTest:

    def setup_class(self, Np=5):
        wrk.clear()
        # Create Topological Network object
        self.net = op.network.Cubic([Np, Np, 1], spacing=1)
        self.net.add_model_collection(
            op.models.collections.geometry.spheres_and_cylinders
        )
        self.net.regenerate_models()
        self.net['pore.diameter'] = 0.5
        self.net['throat.diameter'] = 0.25
        self.phase = op.phase.Air(network=self.net)
        self.inlets = [0]
        self.outlets = [Np*Np - 1]

    def run_mp(self, trapping=False, residual=False, snap=False,
               plot=False, flowrate=None):
        IP_1 = mpc(network=self.net)
        if snap:
            IP_1.settings['snap_off'] = 'throat.snap_off'
        IP_1.setup(phase=self.phase)
        IP_1.set_inlets(pores=self.inlets)
        if residual:
            IP_1.set_residual(pores=self.phase['pore.occupancy'])
        IP_1.run()
        if trapping:
            IP_1.set_outlets(self.outlets)
            IP_1.apply_trapping()
        inv_points = np.arange(0, 100, 1)
        # returns data as well as plotting
        alg_data = IP_1.get_intrusion_data(inv_points=inv_points)
        self.phase.update(IP_1.results(Pc=inv_points.max()))
        if plot:
            plt.figure()
            L = np.sqrt(self.net.Np).astype(int)
            plt.imshow(IP_1['pore.invasion_sequence'].reshape([L, L]),
                       cmap=plt.get_cmap('Blues'))
        if flowrate is not None:
            IP_1.apply_flow(flowrate=flowrate)
        self.alg = IP_1
        return alg_data

    # def test_coop_pore_filling(self):
    #     pn = op.network.Cubic(shape=[3, 3, 3], spacing=2.5e-5)
    #     pn.add_model_collection(op.models.collections.geometry.spheres_and_cylinders)
    #     pn.regenerate_models()
    #     pn['throat.diameter'] = 1.5e-5
    #     pn['pore.diameter'] = 2e-5
    #     pn.add_model(propname='throat.centroid',
    #                  model=op.models.geometry.throat_centroid.pore_coords)
    #     pn.add_model(propname='throat.normal',
    #                  model=op.models.geometry.throat_vector.pore_to_pore)
    #     water = op.phase.Water(network=pn)
    #     water['pore.contact_angle'] = 60
    #     r_tor = 5e-6
    #     water.add_model(propname='throat.entry_pressure',
    #                     model=op.models.physics.meniscus.purcell,
    #                     r_toroid=r_tor,
    #                     mode='max')
    #     water.add_model(propname='throat.meniscus',
    #                     model=op.models.physics.meniscus.purcell,
    #                     mode='men',
    #                     r_toroid=r_tor,
    #                     target_Pc=5000)
    #     water.regenerate_models()
    #     water['pore.entry_pressure'] = 0.0
    #     ip = op.algorithms.MixedInvasionPercolationCoop(network=pn)
    #     ip.setup(phase=water)
    #     ip.setup(cooperative_pore_filling='throat.meniscus')
    #     points = np.arange(0.1, 1, 0.05)*ip._max_pressure()
    #     ip.setup_coop_filling(inv_points=points)
    #     ip.set_inlets(pores=pn.pores('bottom'))
    #     ip.run()
    #     assert np.any(~np.isnan(ip.tt_Pc.data[0]))


if __name__ == '__main__':
    t = MixedPercolationCoopTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
