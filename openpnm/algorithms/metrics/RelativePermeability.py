import numpy as np
from openpnm import models
from openpnm.utils import logging
from openpnm.phases import GenericPhase
from openpnm.algorithms import GenericAlgorithm, StokesFlow
logger = logging.getLogger(__name__)


default_settings = {
    'wp': None,
    'nwp': None,
    'conduit_hydraulic_conductance': 'throat.conduit_hydraulic_conductance',
    'hydraulic_conductance': 'throat.hydraulic_conductance',
    'pore.invasion_sequence': 'pore.invasion_sequence',
    'throat.invasion_sequence': 'throat.invasion_sequence',
    'flow_inlet': None,
    'flow_outlet': None,
    'Snwp_num': None,
}


class RelativePermeability(GenericAlgorithm):
    r"""
    A subclass of Generic Algorithm to calculate relative permeabilities
    of fluids in a drainage process.

    The main roles of this subclass are to get invasion sequence and
    implement a method for calculating the relative permeabilities of the
    fluids flowing in three directions.

    Notes
    -----
    1. The results can be plotted using `plot_Kr_curves`, and numerical
       data can be obtained with `get_Kr_data`.

    2. Properties related to the invading phase have subscript 'nwp',
       while those related to the defending phase (if there is any) are
       named by subscript 'wp'.

    3. As the relative permeability definition is the ratio of effective
       to absolute permeability, domain length and area calculation will
       appear in both nominator and denominator. Ignoring those variables,
       we only use the flow rate of the phase of interest in single and
       multiphase permeability calculation.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update(default_settings)
        self.settings.update(settings)
        self.Kr_values = {'sat': dict(),
                          'relperm_wp': dict(),
                          'relperm_nwp': dict(),
                          'perm_abs_wp': dict(),
                          'perm_abs_nwp': dict(),
                          'results': {'sat': [], 'kr_wp': [], 'kr_nwp': []}}

    def setup(self, invading_phase=None, defending_phase=None,
              invasion_sequence=None, flow_inlets=None, flow_outlets=None,
              Snwp_num=100):
        r"""
        Assigns values to the algorithm ``settings``.

        This part is revised so that it can be applied on either 2D or 3D
        networks. The dimension is found by the number of boundary faces
        existed in the network asumming all boundary faces are labeled.

        Parameters
        ----------
        invading_phase: str
            The invading or non-wetting phase.

        defending_phase: str, optional
            If defending phase is specified, then it's permeability will
            also be calculated, otherwise only the invading phase is
            considered.

        invasion_sequence: str (default is 'invasion_sequence')
            The dictionary key on the invading phase object where the
            invasion sequence is stored. The default from both the IP and
            OP algorithms is 'invasion_sequence', so this is the default
            here.

        flow_inlets: dict, optional
            A dictionary containing the direction of the flow and the
            label of the pores as the inlet pores of that direction.

        flow_outlets: dict, optional
            A dictionary containing the direction of the flow and the
            label of the pores as the outlet pores of that direction.

        swnp_num: int
           The number of saturation points (equidistant points), at which
           the Kr values will be calculated.

        """
        self.settings['Snwp_num'] = Snwp_num
        network = self.project.network
        if invading_phase is not None:
            self.settings['nwp'] = invading_phase
        if defending_phase is not None:
            self.settings['wp'] = defending_phase
        if (invasion_sequence == 'invasion_sequence'):
            nwp = self.project[self.settings['nwp']]
            seq_p = nwp['pore.invasion_sequence']
            self['pore.invasion_sequence'] = seq_p
            seq_t = nwp['throat.invasion_sequence']
            self['throat.invasion_sequence'] = seq_t
            dimension = 0
            x, y, z = False, False, False
            if 'pore.left' in network.keys():
                dimension += 1
                x = True
            if 'pore.front' in network.keys():
                dimension += 1
                y = True
            if 'pore.top' in network.keys():
                dimension += 1
                z = True
            if dimension == 2:
                if (x and y):
                    self.settings['flow_inlets'] = {'x': 'left',
                                                    'y': 'front'}
                    self.settings['flow_outlets'] = {'x': 'right',
                                                     'y': 'back'}
                elif (y and z):
                    self.settings['flow_inlets'] = {'y': 'front',
                                                    'z': 'top'}
                    self.settings['flow_outlets'] = {'y': 'back',
                                                     'z': 'bottom'}
                elif (x and z):
                    self.settings['flow_inlets'] = {'x': 'left',
                                                    'z': 'top'}
                    self.settings['flow_outlets'] = {'x': 'right',
                                                     'z': 'bottom'}
            else:
                self.settings['flow_inlets'] = {'x': 'left',
                                                'y': 'front',
                                                'z': 'top'}
                self.settings['flow_outlets'] = {'x': 'right',
                                                 'y': 'back',
                                                 'z': 'bottom'}
        if flow_inlets is not None:
            for keys in flow_inlets.keys():
                self.settings['flow_inlets'][keys] = flow_inlets[keys]
                self.settings['flow_outlets'][keys] = flow_outlets[keys]

    def _regenerate_models(self):
        r"""
        Updates the multiphase physics model for each saturation.
        """
        prop = self.settings['conduit_hydraulic_conductance']
        prop_q = self.settings['hydraulic_conductance']
        if self.settings['wp'] is not None:
            wp = self.project[self.settings['wp']]
            modelwp = models.physics.multiphase.conduit_conductance
            wp.add_model(model=modelwp, propname=prop,
                         throat_conductance=prop_q, mode='medium')
        nwp = self.project[self.settings['nwp']]
        modelnwp = models.physics.multiphase.conduit_conductance
        nwp.add_model(model=modelnwp, propname=prop,
                      throat_conductance=prop_q, mode='medium')

    def _abs_perm_calc(self, phase, flow_pores):
        r"""
        Calculates absolute permeability of the medium using StokesFlow
        algorithm.

        The direction of flow is defined by flow_pores. This permeability
        is normalized by variables in darcy's law other than the rate.

        Parameters
        ----------
        phase: GenericPhase
            The phase for which the flow rate is calculated.

        flow_pores: np.ndarray
            Boundary pores that will have constant value boundary
            condition to in StokesFlow algorithm. First element is the
            inlet face (pores) for flow of invading phase through porous
            media. Second element is the outlet face (pores).

        Returns
        -------
        K_abs: float
            The value of absolute permeability of the invading phase in
            the direction that is defined by flow_pores.

        Notes
        -----
        Absolute permeability is not dependent to the phase, but here
        we just need to calculate the rate instead of all variables that
        are contributing to the darcy's law.

        """
        network = self.project.network
        St_p = StokesFlow(network=network, phase=phase)
        St_p.set_value_BC(pores=flow_pores[0], values=1)
        St_p.set_value_BC(pores=flow_pores[1], values=0)
        St_p.run()
        val = np.sum(abs(St_p.rate(pores=flow_pores[1])))
        K_abs = val
        self.project.purge_object(obj=St_p)
        return K_abs

    def _eff_perm_calc(self, flow_pores):
        r"""
        Calculates effective permeability of each phase using StokesFlow
        algorithm with updated multiphase physics models to account for
        the multiphase flow.

        The direction of the flow is defined by flow_pores.

        Parameters
        ----------
        flow_pores: np.ndarray
            Boundary pores that will have constant value boundary
            condition in StokesFlow algorithm. First element is the inlet
            face (pores) for flow of invading phase through porous media.
            Second element is the outlet face (pores).

        Returns
        -------
        output: list
            The value of effective permeability of defending (if there is
            any) and invading phase in the direction that is defined by
            flow_pores.

        Notes
        -----
        1. To account for multiphase flow, multiphase physics model is
           added and updated in each saturation (saturation is related to
           the presence of another phase). Here, the conduit_hydraulic
           conductance is used as the conductance required by stokes flow
           algorithm.

        2. All variables except for the rate in darcy's law will be the
           same in relative permeability ratio. The effective rate
           represents the effective permeability in the nominator of
           relative permeability ratio.

        """
        network = self.project.network
        self._regenerate_models()
        if self.settings['wp'] is not None:
            wp = self.project[self.settings['wp']]
            St_mp_wp = StokesFlow(network=network, phase=wp)
            St_mp_wp.setup(conductance='throat.conduit_hydraulic_conductance')
            St_mp_wp.set_value_BC(pores=flow_pores[0], values=1)
            St_mp_wp.set_value_BC(pores=flow_pores[1], values=0)
            St_mp_wp.run()
            Kewp = np.sum(abs(St_mp_wp.rate(pores=flow_pores[1])))
            self.project.purge_object(obj=St_mp_wp)
        else:
            Kewp = None
            pass
        nwp = self.project[self.settings['nwp']]
        St_mp_nwp = StokesFlow(network=network, phase=nwp)
        St_mp_nwp.set_value_BC(pores=flow_pores[0], values=1)
        St_mp_nwp.set_value_BC(pores=flow_pores[1], values=0)
        St_mp_nwp.setup(conductance='throat.conduit_hydraulic_conductance')
        St_mp_nwp.run()
        Kenwp = np.sum(abs(St_mp_nwp.rate(pores=flow_pores[1])))
        Kenwp = Kenwp
        self.project.purge_object(obj=St_mp_nwp)
        return [Kewp, Kenwp]

    def _sat_occ_update(self, i):
        r"""
        Calculates the saturation of each phase using the invasion
        sequence from either invasion percolation or ordinary percolation.

        Parameters
        ----------
        i: int
            The invasion_sequence limit for masking pores/throats that
            have already been invaded within this limit range. The
            saturation is found by adding the volume of pores and thorats
            that meet this sequence limit divided by the bulk volume.

        """
        network = self.project.network
        pore_mask = self['pore.invasion_sequence'] < i
        throat_mask = self['throat.invasion_sequence'] < i
        sat_p = np.sum(network['pore.volume'][pore_mask])
        sat_t = np.sum(network['throat.volume'][throat_mask])
        sat1 = sat_p + sat_t
        bulk = network['pore.volume'].sum() + network['throat.volume'].sum()
        sat = sat1/bulk
        nwp = self.project[self.settings['nwp']]
        nwp['pore.occupancy'] = pore_mask
        nwp['throat.occupancy'] = throat_mask
        if self.settings['wp'] is not None:
            wp = self.project[self.settings['wp']]
            wp['throat.occupancy'] = 1-throat_mask
            wp['pore.occupancy'] = 1-pore_mask
        return sat

    def run(self, Snwp_num=None):
        r"""
        Calculates the saturation of each phase using the invasion
        sequence.

        Parameters
        ----------
        Snwp_num: int, optional
            Number of saturation point to calculate the relative
            permseability values. If not given, default value is 100.
            Saturation points will be Snwp_num (or 100 by default)
            equidistant points in range [0,1].

        Notes
        -----
        1. For three directions of flow the absolute permeability values
           will be calculated using _abs_perm_calc.

        2. For each saturation point, the saturation values are calculated
           by _sat_occ_update. This function also updates occupancies of
           each phase in pores/throats. Effective permeabilities of each
           phase is then calculated. Relative permeability is defined by
           devision of K_eff and K_abs.

        """
        if Snwp_num is None:
            Snwp_num = self.settings['Snwp_num']
        net = self.project.network
        K_dir = set(self.settings['flow_inlets'].keys())
        for dim in K_dir:
            flow_pores = [net.pores(self.settings['flow_inlets'][dim]),
                          net.pores(self.settings['flow_outlets'][dim])]
            if self.settings['wp'] is not None:
                phase = self.project[self.settings['wp']]
                K_abs = self._abs_perm_calc(phase, flow_pores)
                self.Kr_values['perm_abs_wp'].update({dim: K_abs})
                relperm_wp = []
            else:
                relperm_wp = None
            phase = self.project[self.settings['nwp']]
            K_abs = self._abs_perm_calc(phase, flow_pores)
            self.Kr_values['perm_abs_nwp'].update({dim: K_abs})
        for dirs in self.settings['flow_inlets']:
            if self.settings['wp'] is not None:
                relperm_wp = []
            else:
                relperm_wp = None
                pass
            relperm_nwp = []
            max_seq = np.max([np.max(self['pore.invasion_sequence']),
                              np.max(
                              self['throat.invasion_sequence'])])
            start = max_seq//Snwp_num
            stop = max_seq
            step = max_seq//Snwp_num
            Snwparr = []
            flow_pores = [net.pores(self.settings['flow_inlets'][dirs]),
                          net.pores(self.settings['flow_outlets'][dirs])]
            for j in range(start, stop, step):
                sat = self._sat_occ_update(j)
                Snwparr.append(sat)
                [Kewp, Kenwp] = self._eff_perm_calc(flow_pores)
                if self.settings['wp'] is not None:
                    relperm_wp.append(Kewp/self.Kr_values['perm_abs_wp'][dirs])
                relperm_nwp.append(Kenwp/self.Kr_values['perm_abs_nwp'][dirs])
            if self.settings['wp'] is not None:
                self.Kr_values['relperm_wp'].update({dirs: relperm_wp})
            self.Kr_values['relperm_nwp'].update({dirs: relperm_nwp})
            self.Kr_values['sat'].update({dirs: Snwparr})

    def plot_Kr_curves(self, fig=None):  # pragma: no cover
        r"""
        Plots the relative permeability curve of the phase(s) in flow
        direction(s) as Kr vs Saturation points.
        """
        import matplotlib.pyplot as plt

        if fig is None:
            fig, ax = plt.subplots()
        else:
            ax = fig.get_axes()[0]

        for inp in self.settings['flow_inlets']:
            if self.settings['wp'] is not None:
                ax.plot(self.Kr_values['sat'][inp],
                        self.Kr_values['relperm_wp'][inp],
                        'o-', label='Kr_wp'+inp)
                ax.plot(self.Kr_values['sat'][inp],
                        self.Kr_values['relperm_nwp'][inp],
                        '*-', label='Kr_nwp'+inp)
            else:
                ax.plot(self.Kr_values['sat'][inp],
                        self.Kr_values['relperm_nwp'][inp],
                        '*-', label='Kr_nwp'+inp)

        ax.set_xlabel('Snw')
        ax.set_ylabel('Kr')
        ax.set_title('Relative Permability Curves')
        ax.legend()

        return fig

    def get_Kr_data(self):
        r"""
        Returns data points of the relative permeability of the phase(s)
        in flow direction(s) and Saturation points.
        """
        self.Kr_values['results']['sat'] = self.Kr_values['sat']
        if self.settings['wp'] is not None:
            self.Kr_values['results']['kr_wp'] = self.Kr_values['relperm_wp']
        else:
            self.Kr_values['results']['kr_wp'] = None
        self.Kr_values['results']['kr_nwp'] = self.Kr_values['relperm_nwp']
        return self.Kr_values['results']
