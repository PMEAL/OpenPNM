from collections import namedtuple

import numpy as np
from tqdm.auto import tqdm

from openpnm._skgraph.simulations import (
    bond_percolation,
    find_connected_clusters,
    site_percolation,
)
from openpnm.algorithms import Algorithm
from openpnm.utils import Docorator, TypedSet

docstr = Docorator()


__all__ = ['Drainage']


@docstr.get_sections(base='DrainageSettings',
                     sections=['Parameters'])
@docstr.dedent
class DrainageSettings:
    r"""
    Parameters
    ----------
    %(AlgorithmSettings.parameters)s

    throat_entry_pressure : str
        The dictionary key for the pore entry pressure array
    pore_volume : str
        The dictionary key for the pore volume array
    throat_volume : str
        The dictionary key for the throat volume array

    """
    phase = ''
    throat_entry_pressure = 'throat.entry_pressure'
    pore_volume = 'pore.volume'
    throat_volume = 'throat.volume'
    variable_props = TypedSet()


class Drainage(Algorithm):
    """A class to simulate drainage."""

    def __init__(self, phase, name='drainage_?', **kwargs):
        super().__init__(name=name, **kwargs)
        self.settings._update(DrainageSettings())
        self.settings['phase'] = phase.name
        self['pore.bc.inlet'] = False
        self['pore.bc.outlet'] = False
        self.reset()

    def reset(self):
        r"""
        Resets the algorithm's main results so that it can be re-run
        """
        self['pore.invaded'] = False
        self['throat.invaded'] = False
        # self['pore.residual'] = False
        # self['throat.residual'] = False
        self['pore.trapped'] = False
        self['throat.trapped'] = False
        self['pore.invasion_pressure'] = np.inf
        self['throat.invasion_pressure'] = np.inf
        self['pore.invasion_sequence'] = -1
        self['throat.invasion_sequence'] = -1

    def _set_residual(self, pores=None, throats=None, mode='add'):  # pragma: no cover
        raise NotImplementedError("The ability to add residual nwp is not ready yet")
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

    def set_inlet_BC(self, pores, mode='add'):
        r"""
        Specify the pores from which invading fluid will enter the domain

        Parameters
        ----------
        pores : ndarray
            List of pore indices
        mode : str
            How the boundary conditions should be applied. Options are:

            ============ =====================================================
            mode         meaning
            ============ =====================================================
            'add'        (default) Adds the supplied boundary conditions to
                         the given locations. Raises an exception if values
                         of any type already exist in the given locations.
            'overwrite'  Adds supplied boundary conditions to the given
                         locations, including overwriting conditions of the
                         given type or any other type that may be present in
                         the given locations.
            'remove'     Removes boundary conditions of the specified type
                         from the specified locations. If ``bctype`` is not
                         specified then *all* types are removed. If no
                         locations are given then values are remvoed from
                         *all* locations.
            ============ =====================================================
        """
        self.set_BC(pores=pores, bcvalues=True, bctype='inlet', mode=mode)

    def set_outlet_BC(self, pores, mode='add'):
        r"""
        Specify the pores from which defending fluid exits the network.

        This is optional and only required if trapping is to be considered.

        Parameters
        ----------
        pores : ndarray
            The list of outlet pores
        mode : str
            How the boundary conditions should be applied. Options are:

            ============ =====================================================
            mode         meaning
            ============ =====================================================
            'add'        (default) Adds the supplied boundary conditions to
                         the given locations. Raises an exception if values
                         of any type already exist in the given locations.
            'overwrite'  Adds supplied boundary conditions to the given
                         locations, including overwriting conditions of the
                         given type or any other type that may be present in
                         the given locations.
            'remove'     Removes boundary conditions of the specified type
                         from the specified locations. If ``bctype`` is not
                         specified then *all* types are removed. If no
                         locations are given then values are remvoed from
                         *all* locations.
            ============ =====================================================
        """
        self.set_BC(pores=pores, bcvalues=True, bctype='outlet', mode=mode)

    def run(self, pressures=25):
        r"""
        Runs the simulation for the pressure points

        Parameters
        ----------
        pressures : int or ndarray
            The number of pressue steps to apply, or an array of specific
            points

        """
        if isinstance(pressures, int):
            phase = self.project[self.settings.phase]
            hi = 1.25*phase[self.settings.throat_entry_pressure].max()
            low = 0.80*phase[self.settings.throat_entry_pressure].min()
            pressures = np.logspace(np.log10(low), np.log10(hi), pressures)
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
        if np.any(self['pore.bc.outlet']):
            self.apply_trapping()

    def _run_special(self, pressure):
        phase = self.project[self.settings.phase]
        Tinv = phase[self.settings.throat_entry_pressure] <= pressure
        # Tinv += self['throat.invaded']
        # Remove trapped throats from this list, if any
        # Tinv[self['throat.trapped']] = False
        # Perform bond_percolation to label invaded clusters
        s_labels, b_labels = bond_percolation(self.network.conns, Tinv)
        # Remove label from any clusters not connected to the inlets
        s_labels, b_labels = find_connected_clusters(
            b_labels, s_labels, self['pore.bc.inlet'], asmask=False)
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
            clusters = np.unique(s[self['pore.bc.outlet']])
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
        # self['pore.trapped'][self['pore.residual']] = False
        # self['throat.trapped'][self['throat.residual']] = False
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


# %%
# def run_examples():
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    import openpnm as op
    plt.rcParams['figure.facecolor'] = 'darkgrey'
    plt.rcParams['axes.facecolor'] = 'grey'

    np.random.seed(0)
    Nx, Ny, Nz = 10, 10, 1
    pn = op.network.Cubic([Nx, Ny, Nz], spacing=1e-5)
    pn.add_model_collection(op.models.collections.geometry.spheres_and_cylinders)
    pn.regenerate_models()
    nwp = op.phase.Phase(network=pn)
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
    drn.set_inlet_BC(pores=pn.pores('left'))
    pressures = np.logspace(np.log10(0.1e6), np.log10(8e6), 40)
    drn.run(pressures)
    drn.set_outlet_BC(pores=pn.pores('right'))
    drn.apply_trapping()

    # %%
    if 1:
        fig, ax = plt.subplots(1, 1, figsize=[20, 20])
        ax.semilogx(*drn.pc_curve(pressures), 'ro-')
        ax.set_ylim([-.05, 1.05])

    if 1:
        pressures = np.unique(drn['pore.invasion_pressure'])
        n = 6
        p = drn['pore.invasion_pressure'] <= pressures[n]
        t = drn['throat.invasion_pressure'] <= pressures[n]
        ax = op.topotools.plot_coordinates(pn, pores=p, s=500,
                                           color_by=drn['pore.invasion_pressure'])
        ax = op.topotools.plot_connections(pn, throats=t, linewidth=5,
                                           color_by=drn['throat.invasion_pressure'],
                                           ax=ax)
        # t = drn['throat.invasion_pressure'] <= p
        # ax = op.topotools.plot_connections(pn, throats=t, c='k', linestyle='--', ax=ax)

    if 0:
        import openpnm as op
        tseq = drn['throat.invasion_pressure']
        pseq = drn['pore.invasion_pressure']
        Pmax = np.amax(tseq[tseq < np.inf])
        pseq[pseq > Pmax] = Pmax + 1
        tseq[tseq > Pmax] = Pmax + 1
        for j, i in enumerate(tqdm(np.unique(tseq)[:-1])):
            ax = op.topotools.plot_connections(pn, tseq <= i, linestyle='--',
                                               c='r', linewidth=3)
            op.topotools.plot_coordinates(pn, pseq <= i, c='b', marker='x',
                                          markersize=100, ax=ax)
            op.topotools.plot_coordinates(pn, drn['pore.trapped'], c='c',
                                          marker='o', markersize=100, ax=ax)
            plt.savefig(f"{str(j).zfill(3)}.png")
            plt.close()

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
        clusters = np.unique(s[drn['pore.bc.outlet']])
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
