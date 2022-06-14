import numpy as np
from tqdm import tqdm
from collections import namedtuple
from openpnm.core import ModelMixin2
from openpnm.algorithms import Drainage
from openpnm.algorithms._solution import SolutionContainer, PressureScan
from openpnm.utils import Docorator, TypedSet
from openpnm._skgraph.simulations import (
    bond_percolation,
    site_percolation,
    mixed_percolation,
    find_connected_clusters,
    find_trapped_sites,
)


docstr = Docorator()


__all__ = ['Imbibition']


@docstr.get_sections(base='ImbibitionSettings',
                     sections=['Parameters'])
@docstr.dedent
class ImbibitionSettings:
    r"""
    Parameters
    ----------
    %(GenericAlgorithmSettings.parameters)s

    throat_entry_pressure : str
        The dictionary key for the pore entry pressure array
    pore_volume : str
        The dictionary key for the pore volume array
    throat_volume : str
        The dictionary key for the throat volume array

    """
    phase = ''
    quantity = 'pressure'
    throat_entry_pressure = 'throat.entry_pressure'
    pore_entry_pressure = 'pore.entry_pressure'
    pore_volume = 'pore.volume'
    throat_volume = 'throat.volume'
    variable_props = TypedSet()


class Imbibition(Drainage):

    def __init__(self, phase, name='imbibition_#', **kwargs):
        super().__init__(phase=phase, name=name, **kwargs)
        self.settings._update(ImbibitionSettings())
        self.settings['phase'] = phase.name
        self._im = self.project.network.im.tolil()

    def reset(self):
        super().reset()
        # Uninvaded pores/throats are denoted with a -np.inf
        self['pore.invasion_pressure'] = -np.inf
        self['throat.invasion_pressure'] = -np.inf

    def set_residual(self, pores=[], throats=[]):
        super().set_residual(pores=pores, throats=throats)
        self['pore.invasion_pressure'][self['pore.invaded']] = np.inf
        self['throat.invasion_pressure'][self['throat.invaded']] = np.inf

    def run(self, pressures):
        pressures = np.sort(np.array(pressures, ndmin=1))[-1::-1]
        msg = 'Performing imbibition simulation'
        for i, p in enumerate(tqdm(pressures, msg)):
            self._run_special(p)
            pmask = self['pore.invaded'] * (self['pore.invasion_pressure'] == -np.inf)
            self['pore.invasion_pressure'][pmask] = p
            tmask = self['throat.invaded'] * (self['throat.invasion_pressure'] == -np.inf)
            self['throat.invasion_pressure'][tmask] = p
        # If any outlets were specified, evaluate trapping
        if np.any(self['pore.outlets']):
            self.apply_trapping()

    def _run_special(self, pressure):
        phase = self.project[self.settings.phase]
        # Perform Percolation -------------------------------------------------
        Pinv = phase[self.settings.pore_entry_pressure] > pressure
        Tinv = phase[self.settings.throat_entry_pressure] > pressure
        # ax = op.topotools.plot_connections(pn, throats=Tinv)
        # ax = op.topotools.plot_coordinates(pn, pores=Pinv, c='w', s=200, ax=ax)

        # Pre-seed invaded locations with residual, if any
        Pinv += self['pore.invaded']
        Tinv += self['throat.invaded']

        # op.topotools.plot_connections(pn, throats=self['throat.invaded'], c='r', ax=ax)

        # Remove trapped throats from this list, if any
        Pinv[self['pore.trapped']] = False
        Tinv[self['throat.trapped']] = False
        # op.topotools.plot_connections(pn, throats=self['throat.trapped'], c='g', ax=ax)

        # Perform site_percolation to label invaded clusters of pores
        s_labels, b_labels = site_percolation(self.network.conns, Pinv)
        # ax = op.topotools.plot_connections(pn, color_by=(b_labels + 10)*(b_labels >= 0))
        # op.topotools.plot_coordinates(pn, color_by=(s_labels + 10)*(s_labels >= 0), ax=ax, s=200)

        # Remove label from any clusters not connected to the inlets
        s_labels, b_labels = find_connected_clusters(
            b_labels, s_labels, self['pore.inlets'], asmask=False)
        # ax = op.topotools.plot_connections(pn, color_by=(b_labels + 10)*(b_labels >= 0))
        # op.topotools.plot_coordinates(pn, color_by=(s_labels + 10)*(s_labels >= 0), ax=ax, s=200)

        # Mark throats connected to invaded pores as also invaded, if they're small enough
        Pinv = np.where(s_labels >= 0)[0]
        try:
            temp = np.unique(np.hstack(self._im.rows[Pinv]))
            Tinv = Tinv*self.to_mask(throats=temp)
        except ValueError:
            Tinv = []

        # Add result to existing invaded locations
        self['pore.invaded'][Pinv] = True
        self['throat.invaded'][Tinv] = True
        # ax = op.topotools.plot_connections(pn, c='w', linestyle='--')
        # op.topotools.plot_connections(pn, throats=self['throat.invaded'], c='b', ax=ax)
        # op.topotools.plot_coordinates(pn, pores=self['pore.invaded'], c='b', ax=ax, s=200)

    def apply_trapping(self, mode='mixed'):
        pressures = np.unique(self['pore.invasion_pressure'])
        pseq = self['pore.invasion_pressure']
        tseq = self['throat.invasion_pressure']
        for p in pressures:
            if mode == 'orig':
                pmask = self['pore.invasion_pressure'] >= p
                s, b = find_trapped_sites(conns=self.network.conns,
                                          occupied_sites=pmask,
                                          outlets=self['pore.outlets'])
                P_trapped = s >= 0
                T_trapped = P_trapped[self.network.conns].all(axis=1)
                P_trapped[self['pore.residual']] = False
                T_trapped[self['throat.residual']] = False
                # ax = op.topotools.plot_connections(pn, throats=(b >= 0))
                # op.topotools.plot_coordinates(pn, color_by=(s + 10)*(s >= 0), ax=ax, s=200)
                self['pore.trapped'] += P_trapped
                self['throat.trapped'] += T_trapped
            else:
                if mode == 'bond':
                    s, b = bond_percolation(conns=self.network.conns,
                                            occupied_bonds=tseq < p)
                elif mode == 'site':
                    s, b = site_percolation(conns=self.network.conns,
                                            occupied_sites=pseq < p)
                elif mode == 'mixed':
                    s, b = mixed_percolation(conns=self.network.conns,
                                             occupied_sites=pseq < p,
                                             occupied_bonds=tseq < p)
                clusters = np.unique(s[self['pore.outlets']])
                self['pore.trapped'] += np.isin(s, clusters, invert=True)*(s >= 0)
                self['throat.trapped'] += np.isin(b, clusters, invert=True)*(b >= 0)
        self['pore.trapped'][self['pore.residual']] = False
        self['throat.trapped'][self['throat.residual']] = False
        self['pore.invaded'][self['pore.trapped']] = False
        self['throat.invaded'][self['throat.trapped']] = False
        self['pore.invasion_pressure'][self['pore.trapped']] = -np.inf
        self['throat.invasion_pressure'][self['throat.trapped']] = -np.inf

    def pc_curve(self, pressures=None):
        if pressures is None:
            pressures = np.unique(self['pore.invasion_pressure'])
        elif isinstance(pressures, int):
            p = np.unique(self['pore.invasion_pressure'])
            p = p[np.isfinite(p)]
            pressures = np.logspace(np.log10(p.min()/2), np.log10(p.max()*2), pressures)
        else:
            pressures = np.array(pressures)
        pressures = np.sort(pressures)[-1::-1]
        pc = []
        s = []
        Vp = self.network[self.settings.pore_volume]
        Vt = self.network[self.settings.throat_volume]
        for p in pressures:
            Snwp_p = self['pore.invasion_pressure'] >= p
            Snwp_t = self['throat.invasion_pressure'] >= p
            pc.append(p)
            s.append(((Snwp_p*Vp).sum() + (Snwp_t*Vt).sum())/(Vp.sum() + Vt.sum()))
        pc_curve = namedtuple('pc_curve', ('pc', 'snwp'))
        data = pc_curve(np.array(pc), 1-np.array(s))
        return data


if __name__ == "__main__":
    import openpnm as op
    import matplotlib.pyplot as plt
    plt.rcParams['figure.facecolor'] = 'grey'
    plt.rcParams['axes.facecolor'] = 'grey'

    np.random.seed(0)
    Nx, Ny, Nz = 20, 20, 20
    pn = op.network.Cubic([Nx, Ny, Nz], spacing=1e-5)
    pn.add_model_collection(op.models.collections.geometry.spheres_and_cylinders)
    pn.regenerate_models()
    nwp = op.phase.GenericPhase(network=pn)
    nwp['throat.surface_tension'] = 0.480
    nwp['throat.contact_angle'] = 140
    nwp.add_model(propname='throat.entry_pressure',
                  model=op.models.physics.capillary_pressure.washburn)
    nwp.add_model(propname='pore.entry_pressure',
                  model=op.models.physics.capillary_pressure.washburn,
                  contact_angle=140,
                  surface_tension=0.480,
                  diameter='pore.diameter')
    pressures = np.logspace(np.log10(0.1e6), np.log10(8e6), 40)

    # %% Perform Primary Drainage
    drn = op.algorithms.Drainage(network=pn, phase=nwp)
    drn.set_inlets(pores=pn.pores('left'))
    drn.run(pressures)
    drn.set_outlets(pores=pn.pores('right'))
    drn.apply_trapping(mode='site')

    # %% Peform Imbibition
    imb = Imbibition(network=pn, phase=nwp)
    imb.set_inlets(pores=pn.pores('left'))
    imb.set_residual(pores=~drn['pore.invaded'], throats=~drn['throat.invaded'])
    imb.run(pressures)
    imb.set_outlets(pores=pn.pores('right'))
    imb.apply_trapping(mode='orig')

    # %% Perform Secondary Drainage
    drn2 = op.algorithms.Drainage(network=pn, phase=nwp)
    drn2.set_inlets(pores=pn.pores('left'))
    drn2.set_residual(pores=~imb['pore.invaded'], throats=~imb['throat.invaded'])
    drn2.run(pressures)
    drn2.set_outlets(pores=pn.pores('right'))
    drn2.apply_trapping(mode='site')

    # %%
    if 1:
        fig, ax = plt.subplots(1, 1)
        ax.semilogx(*drn.pc_curve(pressures), 'wo-', label='prinmary drainage')
        ax.semilogx(*imb.pc_curve(pressures), 'ko-', label='imbibition')
        ax.semilogx(*drn2.pc_curve(pressures), 'yo-', label='seconary drainage')
        ax.set_ylim([-.05, 1.05])
        ax.set_xlabel('Capillary Pressure [Pa]')
        ax.set_ylabel('Non-wetting Phase Saturation')
        fig.legend()

    # %%
    if 0:
        for p in range(len(pressures)):
            ax = op.topotools.plot_coordinates(
                network=pn, pores=imb['pore.inlets'],
                marker='s', edgecolor='k', c='grey', s=400, label='inlets')
            ax = op.topotools.plot_coordinates(
                network=pn, pores=imb['pore.outlets'],
                ax=ax, marker='d', edgecolor='k', c='grey', s=400, label='outlets')
            ax = op.topotools.plot_connections(
                network=pn, throats=nwp['throat.entry_pressure'] > pressures[p],
                c='white', ax=ax, label='Invadable throats')
            ax = op.topotools.plot_connections(
                network=pn, throats=sol2['throat.invaded'][:, p],
                ax=ax, label='Invaded throats')
            ax = op.topotools.plot_coordinates(
                network=pn, pores=sol2['pore.invaded'][:, p],
                s=100, ax=ax, label='Invaded pores')
            ax = op.topotools.plot_coordinates(
                network=pn, pores=sol2['pore.trapped'][:, p],
                c='green', s=100, ax=ax, label='Trapped pores')
            ax = op.topotools.plot_connections(
                network=pn, throats=sol2['throat.trapped'][:, p],
                c='black', linestyle='--', ax=ax, label='Trapped throats')
            fig = plt.gcf()
            fig.legend(loc='center left', fontsize='large')
            # ax = op.topotools.plot_coordinates(
            #     network=pn, pores=~imb['pore.invaded'],
            #     c='grey', ax=ax)
            plt.savefig(f'figure{p}.png')























