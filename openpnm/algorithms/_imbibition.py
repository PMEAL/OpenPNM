import numpy as np
from openpnm.core import ModelMixin2
from openpnm.algorithms import GenericAlgorithm
from openpnm.algorithms._solution import SolutionContainer, PressureScan
from openpnm.utils import Docorator, TypedSet
from openpnm._skgraph.simulations import (
    bond_percolation,
    site_percolation,
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


class Imbibition(ModelMixin2, GenericAlgorithm):

    def __init__(self, phase, name='imbibition_#', **kwargs):
        super().__init__(name=name, **kwargs)
        self.settings._update(ImbibitionSettings())
        self.settings['phase'] = phase.name
        self._im = self.project.network.im.tolil()
        self.reset()

    def reset(self):
        self['pore.inlets'] = False
        self['pore.outlets'] = False
        self['pore.invaded'] = False
        self['throat.invaded'] = False
        self['pore.residual'] = False
        self['throat.residual'] = False
        self['pore.trapped'] = False
        self['throat.trapped'] = False
        self.soln = SolutionContainer()

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

    def set_residual(self, pores=None, throats=None, mode='add'):
        if pores is not None:
            self['pore.invaded'][pores] = True
            self['pore.residual'][pores] = True
        if throats is not None:
            self['throat.invaded'][throats] = True
            self['throat.residual'][throats] = True

    def set_trapped(self, pores=None, throats=None, mode='add'):
        # A pore/throat cannot be both trapped and invaded so set invaded=False
        if pores is not None:
            self['pore.invaded'][pores] = False
            self['pore.trapped'][pores] = True
        if throats is not None:
            self['throat.invaded'][throats] = False
            self['throat.trapped'][throats] = True

    def run(self, pressures):
        pressures = np.array(pressures, ndmin=1)
        Nx = pressures.size
        self.soln['pore.invaded'] = \
            PressureScan(pressures, np.zeros([self.Np, Nx], dtype=bool))
        self.soln['throat.invaded'] = \
            PressureScan(pressures, np.zeros([self.Nt, Nx], dtype=bool))
        self.soln['pore.trapped'] = \
            PressureScan(pressures, np.zeros([self.Np, Nx], dtype=bool))
        self.soln['throat.trapped'] = \
            PressureScan(pressures, np.zeros([self.Nt, Nx], dtype=bool))
        for item in self.models.keys():
            key = item.split('@')[0]
            element, prop = key.split('.', 1)
            self.soln[key] = \
                PressureScan(pressures, np.zeros([self._count(element), Nx],
                                                 dtype=float))
        for i, p in enumerate(pressures):
            self._run_special(p)
            self.soln['pore.invaded'][:, i] = self['pore.invaded']
            self.soln['throat.invaded'][:, i] = self['throat.invaded']
            self.soln['pore.trapped'][:, i] = self['pore.trapped']
            self.soln['throat.trapped'][:, i] = self['throat.trapped']
            for item in self.models.keys():
                key = item.split('@')[0]
                self.soln[key][:, i] = self[key]
        return self.soln

    def _run_special(self, pressure):
        phase = self.project[self.settings.phase]
        # Perform Percolation -------------------------------------------------
        Pinv = phase[self.settings.pore_entry_pressure] > pressure
        Tinv = phase[self.settings.throat_entry_pressure] > pressure
        ax = op.topotools.plot_connections(pn, throats=Tinv)
        ax = op.topotools.plot_coordinates(pn, pores=Pinv, c='w', s=200, ax=ax)

        # Pre-seed invaded locations with residual, if any
        Pinv += self['pore.invaded']
        Tinv += self['throat.invaded']

        # op.topotools.plot_connections(pn, throats=self['throat.invaded'], c='r', ax=ax)

        # Remove trapped throats from this list, if any
        Pinv[self['pore.trapped']] = False
        Tinv[self['throat.trapped']] = False
        # op.topotools.plot_connections(pn, throats=self['throat.trapped'], c='g', ax=ax)

        # Perform site_percolation to label invaded clusters of pores
        Pinv = self.to_mask(pores=Pinv)
        s_labels, b_labels = site_percolation(self.network.conns, Pinv)
        # ax = op.topotools.plot_connections(pn, color_by=(b_labels + 10)*(b_labels >= 0))
        # op.topotools.plot_coordinates(pn, color_by=(s_labels + 10)*(s_labels >= 0), ax=ax, s=200)

        # Remove label from any clusters not connected to the inlets
        s_labels, b_labels = find_connected_clusters(
            b_labels, s_labels, self['pore.inlets'], asmask=False)
        ax = op.topotools.plot_connections(pn, color_by=(b_labels + 10)*(b_labels >= 0))
        op.topotools.plot_coordinates(pn, color_by=(s_labels + 10)*(s_labels >= 0), ax=ax, s=200)

        # Mark throats connected to invaded pores as also invaded, if they're small enough
        Pinv = np.where(s_labels >= 0)[0]
        Tinv = Tinv*self.to_mask(throats=np.unique(np.hstack(self._im.rows[Pinv])))

        # Add result to existing invaded locations
        self['pore.invaded'][Pinv] = True
        self['throat.invaded'][Tinv] = True
        ax = op.topotools.plot_connections(pn, c='w', linestyle='--')
        op.topotools.plot_connections(pn, throats=self['throat.invaded'], c='b', ax=ax)
        op.topotools.plot_coordinates(pn, pores=self['pore.invaded'], c='b', ax=ax, s=200)

        # Update invasion status and pressure
        pc = self.settings.quantity
        self['pore.'+pc] = pressure*self['pore.invaded']
        self['throat.'+pc] = pressure*self['throat.invaded']
        self.regenerate_models()

        # Trapping ------------------------------------------------------------
        # If any outlets were specified, evaluate trapping
        if np.any(self['pore.outlets']):
            s, b = find_trapped_sites(conns=self.network.conns,
                                      occupied_sites=self['pore.invaded'],
                                      outlets=self['pore.outlets'])
            P_trapped = s >= 0
            T_trapped = P_trapped[self.network.conns].any(axis=1)
            # ax = op.topotools.plot_connections(pn, throats=(b >= 0))
            # op.topotools.plot_coordinates(pn, color_by=(s + 10)*(s >= 0), ax=ax, s=200)
            self['pore.trapped'] += P_trapped
            self['throat.trapped'] += T_trapped

    def pc_curve(self, pressures=None, y='invaded'):
        if pressures is None:
            pressures = self.soln['pore.invaded'].p
        elif isinstance(pressures, int):
            pressures = np.linspace(self.soln['pore.invaded'].p[0],
                                    self.soln['pore.invaded'].p[-1],
                                    pressures)
        else:
            pressures = np.array(pressures)
        pc = []
        s = []
        Snwp_p = self.soln['pore.'+y]
        Snwp_t = self.soln['throat.'+y]
        Vp = self.network[self.settings.pore_volume]
        Vt = self.network[self.settings.throat_volume]
        for p in pressures:
            pc.append(p)
            s.append(((Snwp_p(p)*Vp).sum() + (Snwp_t(p)*Vt).sum())/(Vp.sum() + Vt.sum()))
        return pc, s


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
    plt.rcParams['figure.facecolor'] = 'grey'
    plt.rcParams['axes.facecolor'] = 'grey'

    np.random.seed(0)
    Nx, Ny, Nz = 20, 20, 1
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

    imb = Imbibition(network=pn, phase=nwp)
    imb.set_inlets(pores=pn.pores('left'))
    pressures = np.logspace(np.log10(0.1e6), np.log10(8e6), 40)
    # sol = imb.run(pressures)

    # # %%
    # if 0:
    #     fig, ax = plt.subplots(1, 1)
    #     ax.semilogx(*drn.pc_curve(y='snwp'), 'wo-')
    #     ax.set_ylim([-.05, 1.05])

    # # %%
    # if 0:
    #     p = 22
    #     ax = op.topotools.plot_coordinates(
    #         network=pn, pores=drn['pore.inlets'],
    #         marker='s', edgecolor='k', c='grey', s=400, label='inlets')
    #     ax = op.topotools.plot_coordinates(
    #         network=pn, pores=pn['pore.right'],
    #         ax=ax, marker='d', edgecolor='k', c='grey', s=400, label='outlets')
    #     ax = op.topotools.plot_connections(
    #         network=pn, throats=nwp['throat.entry_pressure'] <= pressures[p],
    #         c='white', ax=ax, label='Invadable throats')
    #     ax = op.topotools.plot_connections(
    #         network=pn, throats=sol['throat.invaded'][:, p],
    #         ax=ax, label='Invaded throats')
    #     ax = op.topotools.plot_coordinates(
    #         network=pn, pores=sol['pore.invaded'][:, p],
    #         s=100, ax=ax, label='Invaded pores')
    #     ax = op.topotools.plot_coordinates(
    #         network=pn, pores=sol['pore.trapped'][:, p],
    #         c='green', s=100, ax=ax, label='Trapped pores')
    #     ax = op.topotools.plot_connections(
    #         network=pn, throats=sol['throat.trapped'][:, p],
    #         c='black', linestyle='--', ax=ax, label='Trapped throats')
    #     fig = plt.gcf()
    #     fig.legend(loc='center left', fontsize='large')
    #     # ax = op.topotools.plot_coordinates(
    #     #     network=pn, pores=~drn['pore.invaded'],
    #     #     c='grey', ax=ax)























