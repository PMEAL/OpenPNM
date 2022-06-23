import numpy as np
from tqdm import tqdm
from collections import namedtuple
from openpnm.algorithms import GenericAlgorithm
from openpnm.utils import Docorator, TypedSet
from openpnm._skgraph.simulations import (
    bond_percolation,
    site_percolation,
    mixed_percolation,
    find_connected_clusters,
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


class Drainage(GenericAlgorithm):

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
        self['pore.invasion_sequence'] = -1
        self['throat.invasion_sequence'] = -1

    def set_residual(self, pores=None, throats=None, mode='add'):
        if pores is not None:
            self['pore.invaded'][pores] = True
            self['pore.residual'][pores] = True
            self['pore.invasion_pressure'][self['pore.invaded']] = -np.inf
            self['pore.invasion_sequence'][pores] = 0
        if throats is not None:
            self['throat.invaded'][throats] = True
            self['throat.residual'][throats] = True
            self['throat.invasion_pressure'][self['throat.invaded']] = -np.inf
            self['throat.invasion_sequence'][throats] = 0

    def set_inlets(self, pores, mode='add'):
        if mode == 'add':
            self['pore.inlets'][pores] = True
        elif mode == 'remove':
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
        elif mode == 'remove':
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
            self['pore.invasion_sequence'][pmask] = i
            tmask = self['throat.invaded'] * (self['throat.invasion_pressure'] == np.inf)
            self['throat.invasion_pressure'][tmask] = p
            self['throat.invasion_sequence'][tmask] = i
        # If any outlets were specified, evaluate trapping
        if np.any(self['pore.outlets']):
            self.apply_trapping()

    def _run_special(self, pressure):
        phase = self.project[self.settings.phase]
        Tinv = phase[self.settings.throat_entry_pressure] <= pressure
        Tinv += self['throat.invaded']
        # Remove trapped throats from this list, if any
        Tinv[self['throat.trapped']] = False
        # Perform bond_percolation to label invaded clusters
        s_labels, b_labels = bond_percolation(self.network.conns, Tinv)
        # Remove label from any clusters not connected to the inlets
        s_labels, b_labels = find_connected_clusters(
            b_labels, s_labels, self['pore.inlets'], asmask=False)
        # Add result to existing invaded locations
        self['pore.invaded'][s_labels >= 0] = True
        self['throat.invaded'][b_labels >= 0] = True

    def apply_trapping(self):
        r"""
        Adjusts the invasion history of pores and throats that are trapped.

        Returns
        -------
        This function returns nothing, but the following adjustments are made
        to the data on the object for trapped pores and throats:

        * ``'pore/throat.trapped'`` is set to ``True``
        * ``'pore/throat.invaded'`` is set to ``False``
        * ``'pore/throat.invasion_pressure'`` is set to ``np.inf``
        * ``'pore/throat.invasion_sequence'`` is set to ``0``

        Notes
        -----
        This search proceeds by the following 3 steps:

        1. A site percolation is applied to *uninvaded* pores and they are set
        to trapped if they belong to a cluster that is not connected to the
        outlets.

        2. All throats which were invaded at a pressure *higher* than either
        of its two neighboring pores are set to trapped, regardless of
        whether the pores themselves are trapped.

        3. All throats which are connected to trapped pores are set to trapped
        as these cannot be invaded since the fluid they contain cannot escape.

        """
        pseq = self['pore.invasion_pressure']
        tseq = self['throat.invasion_pressure']
        # Firstly, find any throats who were invaded at a pressure higher than
        # either of its two neighboring pores
        temp = (pseq[self.network.conns].T > tseq).T
        self['throat.trapped'][np.all(temp, axis=1)] = True
        # Now scan through and use site percolation to find other trapped
        # clusters of pores
        for p in np.unique(pseq):
            s, b = site_percolation(conns=self.network.conns,
                                    occupied_sites=pseq > p)
            # Identify cluster numbers connected to the outlets
            clusters = np.unique(s[self['pore.outlets']])
            # Find ALL throats connected to any trapped site, since these
            # throats must also be trapped, and update their cluster numbers
            Ts = self.network.find_neighbor_throats(pores=s >= 0)
            b[Ts] = np.amax(s[self.network.conns], axis=1)[Ts]
            # Finally, mark pores and throats as trapped if their cluster
            # numbers are NOT connected to the outlets
            self['pore.trapped'] += np.isin(s, clusters, invert=True)*(s >= 0)
            self['throat.trapped'] += np.isin(b, clusters, invert=True)*(b >= 0)
        # Use the identified trapped pores and throats to update the other
        # data on the object accordingly
        self['pore.trapped'][self['pore.residual']] = False
        self['throat.trapped'][self['throat.residual']] = False
        self['pore.invaded'][self['pore.trapped']] = False
        self['throat.invaded'][self['throat.trapped']] = False
        self['pore.invasion_pressure'][self['pore.trapped']] = np.inf
        self['throat.invasion_pressure'][self['throat.trapped']] = np.inf
        self['pore.invasion_sequence'][self['pore.trapped']] = -1
        self['throat.invasion_sequence'][self['throat.trapped']] = -1

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
    Nx, Ny, Nz = 50, 50, 1
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
    # drn.set_residual(pores=np.random.randint(0, Nx*Ny, int(Nx*Ny/10)))
    drn.set_inlets(pores=pn.pores('left'))
    pressures = np.logspace(np.log10(0.1e6), np.log10(8e6), 40)
    drn.run(pressures)
    drn.set_outlets(pores=pn.pores('right'))
    drn.apply_trapping()

    # %%
    if 0:
        fig, ax = plt.subplots(1, 1, figsize=[20, 20])
        ax.semilogx(*drn.pc_curve(pressures), 'ro-')
        ax.set_ylim([-.05, 1.05])

    if 0:
        import openpnm as op
        tseq = drn['throat.invasion_pressure']
        pseq = drn['pore.invasion_pressure']
        Pmax = np.amax(tseq[tseq < np.inf])
        pseq[pseq > Pmax] = Pmax + 1
        tseq[tseq > Pmax] = Pmax + 1
        for j, i in enumerate(tqdm(np.unique(tseq)[:-1])):
            ax = op.topotools.plot_connections(pn, tseq<=i, linestyle='--', c='r', linewidth=3)
            op.topotools.plot_coordinates(pn, pseq<=i, c='b', marker='x', markersize=100, ax=ax)
            op.topotools.plot_coordinates(pn, drn['pore.trapped'], c='c', marker='o', markersize=100, ax=ax)
            plt.savefig(f"{str(j).zfill(3)}.png")
            plt.close()

    # %%
    if 0:
        drn.reset()
        p = .9e6
        ps = np.logspace(np.log10(1e5), np.log10(p), 20)
        drn.run(ps)
        # drn.apply_trapping()
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

    # %%
    if 0:
        drn = Drainage(network=pn, phase=nwp)
        drn.set_inlets(pores=pn.pores('left'))
        pressures = np.logspace(np.log10(0.1e6), np.log10(8e6), 40)
        drn.run(pressures)
        drn.set_outlets(pores=pn.pores('right'))
        pressures = np.unique(drn['pore.invasion_pressure'])
        pseq = drn['pore.invasion_pressure']
        p = pressures[4]
        s, b = site_percolation(conns=pn.conns, occupied_sites=pseq > p)
        clusters = np.unique(s[drn['pore.outlets']])
        Ts = pn.find_neighbor_throats(pores=s >= 0)
        b[Ts] = np.amax(s[pn.conns], axis=1)[Ts]
        drn['pore.trapped'] += np.isin(s, clusters, invert=True)*(s >= 0)
        drn['throat.trapped'] += np.isin(b, clusters, invert=True)*(b >= 0)
        ax = op.topotools.plot_coordinates(pn, pores=drn['pore.trapped'],
                                           color_by=s)
        ax = op.topotools.plot_coordinates(pn, pores=pseq <= p, c='k', ax=ax)
        ax = op.topotools.plot_connections(pn, throats=drn['throat.trapped'],
                                           color_by=b, ax=ax)
        t = drn['throat.invasion_pressure'] <= p
        ax = op.topotools.plot_connections(pn, throats=t, c='k', linestyle='--', ax=ax)
