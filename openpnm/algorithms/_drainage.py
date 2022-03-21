import numpy as np
from openpnm.algorithms import GenericAlgorithm
from openpnm._skgraph.simulations import bond_percolation, trim_disconnected_clusters
from openpnm.topotools import remove_isolated_clusters, ispercolating
from openpnm.utils import logging, SettingsAttr, Docorator


__all__ = ['Drainage']


docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sections(base='DrainageSettings',
                     sections=['Parameters'])
@docstr.dedent
class DrainageSettings:
    r"""
    Parameters
    ----------
    %(GenericAlgorithmSettings.parameters)s
    access_limited : bool
        If ``True`` then invading fluid must be connected to the specified
        inlets
    mode : str
        Controls whether pore or throat entry threshold values are used.
        Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'site'       the pore entry is considered
            'bond'       the throat values are considered.
            ===========  =====================================================

    pore_entry_threshold : str
        The dictionary key for the pore entry pressure array
    throat_entry_threshold : str
        The dictionary key for the pore entry pressure array
    pore_volume : str
        The dictionary key for the pore volume array
    throat_volume : str
        The dictionary key for the throat volume array

    """
    phase = ''
    throat_entry_threshold = 'throat.entry_pressure'
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


    def set_inlets(self, pores, mode='add'):
        r"""
        Applies invading phase pressure boundary conditions to the specified pores

        Parameters
        ----------
        pores : array_like
            The pores where the given pressure should be applied.
        mode : str
            Controls how the BCs are added. Options are:

            ============ ==========================================================
            mode         description
            ============ ==========================================================
            'add'        Adds the given locations to any already present
            ============ ==========================================================

        """
        self['pore.inlets'][pores] = True


    def run(self, pressure):
        r"""
        Applies the given pressure to the existing inlets
        """
        phase = self.project[self.settings.phase]
        Tinv = phase[self.settings.throat_entry_threshold] <= pressure
        s_labels, b_labels = bond_percolation(self.network.conns, Tinv)
        s_labels, b_labels = trim_disconnected_clusters(b_labels, s_labels, self['pore.inlets'])
        self['pore.invaded']  = s_labels
        self['throat.invaded']  = b_labels


if __name__ == "__main__":
    import openpnm as op
    pn = op.network.Cubic([20, 20, 1], spacing=1e-5)
    geo = op.geometry.SpheresAndCylinders(network=pn, pores=pn.Ps, throats=pn.Ts)
    nwp = op.phase.GenericPhase(network=pn)
    nwp['throat.surface_tension'] = 0.480
    nwp['throat.contact_angle'] = 180
    phys = op.physics.GenericPhysics(network=pn, phase=nwp, geometry=geo)
    phys.add_model(propname='throat.entry_pressure',
                   model=op.models.physics.capillary_pressure.washburn)

    P = .8e6
    drn = Drainage(network=pn, phase=nwp)
    drn.set_inlets(pores=pn.pores('left'))
    drn.run(pressure=P)

    ax = op.topotools.plot_coordinates(pn, pores=drn['pore.inlets'])
    ax = op.topotools.plot_coordinates(pn, pores=pn.pores('right'), c='grey', ax=ax)
    ax = op.topotools.plot_connections(pn, throats=nwp['throat.entry_pressure'] <= P, c='grey', ax=ax)
    ax = op.topotools.plot_connections(pn, throats=drn['throat.invaded'], ax=ax)
    ax = op.topotools.plot_coordinates(pn, pores=drn['pore.invaded'], ax=ax)
    ax = op.topotools.plot_coordinates(pn, pores=~drn['pore.invaded'], c='grey', ax=ax)































