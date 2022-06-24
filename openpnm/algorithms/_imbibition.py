import numpy as np
from tqdm import tqdm
from collections import namedtuple
from openpnm.algorithms import Drainage
from openpnm.utils import Docorator, TypedSet
from openpnm._skgraph.simulations import (
    site_percolation,
    find_connected_clusters,
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

    def _set_residual(self, pores=[], throats=[]):
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
        Pinv = phase[self.settings.pore_entry_pressure] > pressure
        Tinv = phase[self.settings.throat_entry_pressure] > pressure

        # Pre-seed invaded locations with residual, if any
        Pinv += self['pore.invaded']
        Tinv += self['throat.invaded']

        # Remove trapped throats from this list, if any
        Pinv[self['pore.trapped']] = False
        Tinv[self['throat.trapped']] = False

        # Perform site_percolation to label invaded clusters of pores
        s_labels, b_labels = site_percolation(self.network.conns, Pinv)

        # Remove label from any clusters not connected to the inlets
        s_labels, b_labels = find_connected_clusters(
            b_labels, s_labels, self['pore.inlets'], asmask=False)

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

    def apply_trapping(self):
        pseq = self['pore.invasion_pressure']
        tseq = self['throat.invasion_pressure']
        for p in pseq:
            s, b = site_percolation(conns=self.network.conns,
                                    occupied_sites=pseq < p)
            clusters = np.unique(s[self['pore.outlets']])
            # Ts = self.network.find_neighbor_throats(pores=s >= 0)
            # b[Ts] = np.amax(s[self.network.conns], axis=1)[Ts]
            self['pore.trapped'] += np.isin(s, clusters, invert=True)*(s >= 0)
            self['throat.trapped'] += np.isin(b, clusters, invert=True)*(b >= 0)
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
    pressures = np.logspace(np.log10(0.1e6), np.log10(2e7), 40)

    # %% Perform Primary Drainage
    drn = op.algorithms.Drainage(network=pn, phase=nwp)
    drn.set_inlets(pores=pn.pores('left'))
    drn.run(pressures)
    drn.set_outlets(pores=pn.pores('right'))
    drn.apply_trapping()

    # %% Peform Imbibition
    imb = Imbibition(network=pn, phase=nwp)
    imb.set_inlets(pores=pn.pores('right'))
    imb.run(pressures)
    imb.set_outlets(pores=pn.pores('left'))
    imb.apply_trapping()

    # %%
    if 1:
        fig, ax = plt.subplots(1, 1)
        ax.semilogx(*drn.pc_curve(pressures), 'wo-', label='prinmary drainage')
        ax.semilogx(*imb.pc_curve(pressures), 'ko-', label='imbibition')
        ax.set_ylim([-.05, 1.05])
        ax.set_xlabel('Capillary Pressure [Pa]')
        ax.set_ylabel('Non-wetting Phase Saturation')
        fig.legend()

    if 0:
        pressures = np.unique(imb['pore.invasion_pressure'])
        n = 4
        p = imb['pore.invasion_pressure'] >= pressures[n]
        t = imb['throat.invasion_pressure'] >= pressures[n]
        ax = op.topotools.plot_coordinates(pn, pores=p, s=500,
                                           color_by=imb['pore.invasion_pressure'])
        ax = op.topotools.plot_connections(pn, throats=t, linewidth=5,
                                           color_by=imb['throat.invasion_pressure'], ax=ax)























