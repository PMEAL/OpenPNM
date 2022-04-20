import numpy as np
from openpnm.algorithms import GenericAlgorithm
from openpnm.algorithms._solution import SolutionContainer, PressureScan
from openpnm._skgraph.simulations import bond_percolation
from openpnm._skgraph.simulations import find_connected_clusters
from openpnm._skgraph.simulations import find_trapped_clusters
from openpnm.topotools import remove_isolated_clusters, ispercolating
from openpnm.utils import SettingsAttr, Docorator


__all__ = ['Drainage']
docstr = Docorator()


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
    throat_entry_pressure = 'throat.entry_pressure'
    pore_volume = 'pore.volume'
    throat_volume = 'throat.volume'


class Drainage(GenericAlgorithm):

    def __init__(self, phase, settings=None, **kwargs):
        self.settings = SettingsAttr(DrainageSettings, settings)
        super().__init__(settings=self.settings, **kwargs)
        self.settings['phase'] = phase.name
        self.reset()

    def reset(self):
        self['pore.inlets'] = False
        self['pore.outlets'] = False
        self['pore.invaded'] = False
        self['throat.invaded'] = False
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
        if throats is not None:
            self['throat.invaded'][throats] = True

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
            PressureScan(pressures, np.zeros([self.Np, Nx], dtype=int))
        self.soln['throat.invaded'] = \
            PressureScan(pressures, np.zeros([self.Nt, Nx], dtype=int))
        self.soln['pore.trapped'] = \
            PressureScan(pressures, np.zeros([self.Np, Nx], dtype=int))
        self.soln['throat.trapped'] = \
            PressureScan(pressures, np.zeros([self.Nt, Nx], dtype=int))
        for i, p in enumerate(pressures):
            self._run_special(p)
            self.soln['pore.invaded'][:, i] = self['pore.invaded']
            self.soln['throat.invaded'][:, i] = self['throat.invaded']
            self.soln['pore.trapped'][:, i] = self['pore.trapped']
            self.soln['throat.trapped'][:, i] = self['throat.trapped']
        return self.soln

    def _run_special(self, pressure):
        phase = self.project[self.settings.phase]
        Tinv = phase[self.settings.throat_entry_pressure] <= pressure
        # Update invaded locations with any residual
        Tinv += self['throat.invaded']
        # Removed trapped throats from this list, if any
        Tinv[self['throat.trapped']] = False
        # Performed bond_percolation to label invaded clusters
        s_labels, b_labels = bond_percolation(self.network.conns, Tinv)
        # Remove label from clusters not connected to the inlets
        s_labels, b_labels = find_connected_clusters(
            b_labels, s_labels, self['pore.inlets'], asmask=False)
        # Add result to existing invaded locations
        self['pore.invaded'][s_labels >= 0] = True
        self['throat.invaded'][b_labels >= 0] = True
        # If any outlets were specified, evaluate trapping
        if np.any(self['pore.outlets']):
            s, b = find_trapped_clusters(conns=self.network.conns,
                                         occupied_bonds=self['throat.invaded'],
                                         outlets=self['pore.outlets'])
            self['pore.trapped'][s >= 0] = True
            self['throat.trapped'][b >= 0] = True

    def pc_curve(self, pressures=None):
        if pressures is None:
            pressures = self.soln['pore.invaded'].p
        elif isinstance(pressures, int):
            pressures = np.linspace(self.soln['pore.invaded'].p[0],
                                    self.soln['pore.invaded'].p[-1],
                                    pressures)
        else:
            pressures = np.array(pressures)
        pc = []
        snwp = []
        Vp = self.soln['pore.invaded']
        Vt = self.soln['throat.invaded']
        for p in pressures:
            pc.append(p)
            snwp.append((Vp(p).sum() + Vt(p).sum())/(Vp(p).size + Vt(p).size))
        return pc, snwp


if __name__ == "__main__":
    import openpnm as op
    import matplotlib.pyplot as plt
    plt.rcParams['figure.facecolor'] = 'grey'
    plt.rcParams['axes.facecolor'] = 'grey'
    np.random.seed(0)
    pn = op.network.Cubic([20, 20, 1], spacing=1e-5)
    geo = op.geometry.SpheresAndCylinders(network=pn, pores=pn.Ps, throats=pn.Ts)
    nwp = op.phase.GenericPhase(network=pn)
    nwp['throat.surface_tension'] = 0.480
    nwp['throat.contact_angle'] = 140
    phys = op.physics.GenericPhysics(network=pn, phase=nwp, geometry=geo)
    phys.add_model(propname='throat.entry_pressure',
                   model=op.models.physics.capillary_pressure.washburn)

    drn = Drainage(network=pn, phase=nwp)
    drn.set_residual(pores=np.random.randint(0, 399, 150))
    drn.set_inlets(pores=pn.pores('left'))
    drn.set_outlets(pores=pn.pores('right'))
    pressures = np.linspace(0.4e6, 2e6, 20)
    sol = drn.run(pressures)

    # %%
    ax = op.topotools.plot_coordinates(pn, pores=drn['pore.inlets'], c='m')
    ax = op.topotools.plot_coordinates(pn, pores=pn['pore.right'], ax=ax, c='m')
    ax = op.topotools.plot_connections(pn, throats=nwp['throat.entry_pressure'] <= pressures[4], c='white', ax=ax)
    ax = op.topotools.plot_connections(pn, throats=drn['throat.invaded'], ax=ax)
    ax = op.topotools.plot_coordinates(pn, pores=drn['pore.invaded'], ax=ax)
    ax = op.topotools.plot_coordinates(pn, pores=drn['pore.trapped'], c='green', ax=ax)
    # ax = op.topotools.plot_coordinates(pn, pores=~drn['pore.invaded'], c='grey', ax=ax)

    # %%
    # ax = op.topotools.plot_connections(pn, throats=drn['throat.invaded'], c='grey', ax=ax)
    # ax = op.topotools.plot_coordinates(pn, pores=drn['pore.invaded'], c='blue', ax=ax)
