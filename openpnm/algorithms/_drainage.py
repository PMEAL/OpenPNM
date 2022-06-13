import numpy as np
from tqdm import tqdm
from collections import namedtuple
from openpnm.core import ModelMixin2
from openpnm.algorithms import GenericAlgorithm
from openpnm.utils import Docorator, TypedSet
from openpnm._skgraph.simulations import (
    bond_percolation,
    site_percolation,
    mixed_percolation,
    find_connected_clusters,
    find_trapped_sites,
    find_trapped_bonds,
)


docstr = Docorator()


__all__ = ['Drainage']


@docstr.get_sections(base='DrainageSettings',
                     sections=['Parameters'])
@docstr.dedent
class DrainageSettings:
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
    pore_volume = 'pore.volume'
    throat_volume = 'throat.volume'
    variable_props = TypedSet()


class Drainage(ModelMixin2, GenericAlgorithm):

    def __init__(self, phase, name='drainage_#', **kwargs):
        super().__init__(name=name, **kwargs)
        self.settings._update(DrainageSettings())
        self.settings['phase'] = phase.name
        self['pore.inlets'] = False
        self['pore.outlets'] = False
        self.reset()

    def reset(self):
        self['pore.invaded'] = False
        self['throat.invaded'] = False
        self['pore.residual'] = False
        self['throat.residual'] = False
        self['pore.trapped'] = False
        self['throat.trapped'] = False
        self['pore.invasion_pressure'] = np.inf
        self['throat.invasion_pressure'] = np.inf

    def set_residual(self, pores=None, throats=None, mode='add'):
        if pores is not None:
            self['pore.invaded'][pores] = True
            self['pore.residual'][pores] = True
            self['pore.invasion_pressure'][self['pore.invaded']] = -np.inf
        if throats is not None:
            self['throat.invaded'][throats] = True
            self['throat.residual'][throats] = True
            self['throat.invasion_pressure'][self['throat.invaded']] = -np.inf

    def set_inlets(self, pores, mode='add'):
        if mode == 'add':
            self['pore.inlets'][pores] = True
        elif mode == 'drop':
            self['pore.inlets'][pores] = False
        elif mode == 'clear':
            self['pore.inlets'] = False
        elif mode == 'overwrite':
            self['pore.inlets'] = False
            self['pore.inlets'][pores] = True
        else:
            raise Exception(f'Unrecognized mode {mode}')

    def set_outlets(self, pores, mode='add'):
        if mode == 'add':
            self['pore.outlets'][pores] = True
        elif mode == 'drop':
            self['pore.outlets'][pores] = False
        elif mode == 'clear':
            self['pore.outlets'] = False
        elif mode == 'overwrite':
            self['pore.outlets'] = False
            self['pore.outlets'][pores] = True
        else:
            raise Exception(f'Unrecognized mode {mode}')

    def run(self, pressures):
        pressures = np.array(pressures, ndmin=1)
        msg = 'Performing drainage simulation'
        for i, p in enumerate(tqdm(pressures, msg)):
            self._run_special(p)
            pmask = self['pore.invaded'] * (self['pore.invasion_pressure'] == np.inf)
            self['pore.invasion_pressure'][pmask] = p
            tmask = self['throat.invaded'] * (self['throat.invasion_pressure'] == np.inf)
            self['throat.invasion_pressure'][tmask] = p
        # If any outlets were specified, evaluate trapping
        if np.any(self['pore.outlets']):
            self.apply_trapping()

    def _run_special(self, pressure):
        phase = self.project[self.settings.phase]
        # Perform Percolation -------------------------------------------------
        Tinv = phase[self.settings.throat_entry_pressure] <= pressure
        # ax = op.topotools.plot_connections(pn, throats=Tinv)

        # Pre-seed invaded locations with residual, if any
        Tinv += self['throat.invaded']
        # op.topotools.plot_connections(pn, throats=self['throat.invaded'], c='r', ax=ax)

        # Remove trapped throats from this list, if any
        Tinv[self['throat.trapped']] = False
        # op.topotools.plot_connections(pn, throats=self['throat.trapped'], c='g', ax=ax)

        # Perform bond_percolation to label invaded clusters
        s_labels, b_labels = bond_percolation(self.network.conns, Tinv)
        # ax = op.topotools.plot_connections(pn, color_by=(b_labels + 10)*(b_labels >= 0))
        # op.topotools.plot_coordinates(pn, color_by=(s_labels + 10)*(s_labels >= 0), ax=ax, s=200)

        # Remove label from any clusters not connected to the inlets
        s_labels, b_labels = find_connected_clusters(
            b_labels, s_labels, self['pore.inlets'], asmask=False)
        # ax = op.topotools.plot_connections(pn, color_by=(b_labels + 10)*(b_labels >= 0))
        # op.topotools.plot_coordinates(pn, color_by=(s_labels + 10)*(s_labels >= 0), ax=ax, s=200)

        # Add result to existing invaded locations
        self['pore.invaded'][s_labels >= 0] = True
        self['throat.invaded'][b_labels >= 0] = True
        # ax = op.topotools.plot_connections(pn, c='w', linestyle='--')
        # op.topotools.plot_connections(pn, throats=self['throat.invaded'], c='b', ax=ax)
        # op.topotools.plot_coordinates(pn, pores=self['pore.invaded'], c='b', ax=ax, s=200)

    def apply_trapping(self, mode='mixed'):
        pressures = np.unique(self['pore.invasion_pressure'])
        tseq = self['throat.invasion_pressure']
        pseq = self['pore.invasion_pressure']
        for p in pressures:
            if mode == 'bond':
                s, b = bond_percolation(conns=self.network.conns,
                                        occupied_bonds=tseq > p)
            elif mode == 'site':
                s, b = site_percolation(conns=self.network.conns,
                                        occupied_sites=pseq > p)
            elif mode == 'mixed':
                s, b = mixed_percolation(conns=self.network.conns,
                                         occupied_sites=pseq > p,
                                         occupied_bonds=tseq > p)
            clusters = np.unique(s[self['pore.outlets']])
            self['pore.trapped'] += np.isin(s, clusters, invert=True)*(s >= 0)
            self['throat.trapped'] += np.isin(b, clusters, invert=True)*(b >= 0)
        self['pore.trapped'][self['pore.residual']] = False
        self['throat.trapped'][self['throat.residual']] = False
        self['pore.invaded'][self['pore.trapped']] = False
        self['throat.invaded'][self['throat.trapped']] = False
        self['pore.invasion_pressure'][self['pore.trapped']] = np.inf
        self['throat.invasion_pressure'][self['throat.trapped']] = np.inf

    def pc_curve(self, pressures=None):
        if pressures is None:
            pressures = np.unique(self['pore.invasion_pressure'])
        elif isinstance(pressures, int):
            p = np.unique(self['pore.invasion_pressure'])
            p = p[np.isfinite(p)]
            pressures = np.logspace(np.log10(p.min()/2), np.log10(p.max()*2), pressures)
        else:
            pressures = np.array(pressures)
        pc = []
        s = []
        Vp = self.network[self.settings.pore_volume]
        Vt = self.network[self.settings.throat_volume]
        for p in pressures:
            Snwp_p = self['pore.invasion_pressure'] <= p
            Snwp_t = self['throat.invasion_pressure'] <= p
            pc.append(p)
            s.append(((Snwp_p*Vp).sum() + (Snwp_t*Vt).sum())/(Vp.sum() + Vt.sum()))
        pc_curve = namedtuple('pc_curve', ('pc', 'snwp'))
        data = pc_curve(np.array(pc), np.array(s))
        return data


def late_filling(target,
                 pc='pore.pressure',
                 pc_star='pore.pc_star',
                 eta=1,
                 swp_star=0.2):
    r"""
    Computes the saturation of each pore at the given pressure

    Parameters
    ----------
    target : dict
        The algorithm dictionary
    pnwp : str
        The name of the array containing the capillary pressure defined as the
        pressure difference between the non-wetting and wetting phases. A value
        less than 0 indicates that the element is not invaded.

    """
    pc = target[pc]
    pc_star = target[pc_star]
    swp = np.ones_like(pc)
    # Skip calc for
    mask = pc > pc_star
    if np.any(mask):
        temp = swp_star*(pc_star/pc)**eta
        swp[mask] = temp[mask]
    mask = (pc < pc_star) * (pc > 0)
    swp[mask] = swp_star
    snwp = 1 - swp
    return snwp


if __name__ == "__main__":
    import openpnm as op
    import matplotlib.pyplot as plt
    plt.rcParams['figure.facecolor'] = 'darkgrey'
    plt.rcParams['axes.facecolor'] = 'grey'

    np.random.seed(0)
    Nx, Ny, Nz = 60, 60, 1
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
    # %%
    drn = Drainage(network=pn, phase=nwp)
    drn.set_residual(pores=np.random.randint(0, Nx*Ny, int(Nx*Ny/10)))
    drn.set_inlets(pores=pn.pores('left'))
    pressures = np.logspace(np.log10(0.1e6), np.log10(8e6), 40)
    drn.run(pressures)
    drn.set_outlets(pores=pn.pores('right'))
    drn.apply_trapping(mode='bond')

    # %%
    if 1:
        fig, ax = plt.subplots(1, 1)
        ax.semilogx(*drn.pc_curve(pressures), 'wo-')
        ax.set_ylim([-.05, 1.05])

    # %%
    if 0:
        drn.reset()
        p = .9e6
        ps = np.logspace(np.log10(1e5), np.log10(p), 20)
        # drn.run(ps)
        # drn.apply_trapping()
        for p in ps:
            drn.run(p)
            drn.apply_trapping()
        ax = op.topotools.plot_coordinates(
            network=pn, pores=drn['pore.inlets'],
            marker='s', edgecolor='k', c='grey', s=400, label='inlets')
        ax = op.topotools.plot_coordinates(
            network=pn, pores=pn['pore.right'],
            ax=ax, marker='d', edgecolor='k', c='grey', s=400, label='outlets')
        ax = op.topotools.plot_connections(
            network=pn, throats=nwp['throat.entry_pressure'] <= p,
            c='white', ax=ax, label='Invadable throats')
        ax = op.topotools.plot_connections(
            network=pn, throats=drn['throat.invaded'],
            ax=ax, label='Invaded throats')
        ax = op.topotools.plot_coordinates(
            network=pn, pores=drn['pore.invaded'],
            s=100, ax=ax, label='Invaded pores')
        ax = op.topotools.plot_coordinates(
            network=pn, pores=drn['pore.trapped'],
            c='green', s=100, ax=ax, label='Trapped pores')
        ax = op.topotools.plot_connections(
            network=pn, throats=drn['throat.trapped'],
            c='black', linestyle='--', ax=ax, label='Trapped throats')
        fig = plt.gcf()
        fig.legend(loc='center left', fontsize='large')
        # ax = op.topotools.plot_coordinates(
        #     network=pn, pores=~drn['pore.invaded'],
        #     c='grey', ax=ax)























