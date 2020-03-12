import numpy as np
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class AdvectionDiffusion(ReactiveTransport):
    r"""
    A subclass of GenericTransport to simulate advection-diffusion

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'quantity': 'pore.concentration',
                   'conductance': 'throat.ad_dif_conductance',
                   'diffusive_conductance': 'throat.diffusive_conductance',
                   'hydraulic_conductance': 'throat.hydraulic_conductance',
                   'pressure': 'pore.pressure',
                   's_scheme': 'exponential',
                   'gui': {'setup':        {'phase': None,
                                            'quantity': '',
                                            'conductance': ''},
                           'set_rate_BC':  {'pores': None,
                                            'values': None},
                           'set_value_BC': {'pores': None,
                                            'values': None},
                           'set_source':   {'pores': None,
                                            'propname': ''}
                           }
                   }
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)
        # Make "conductance" iterative_prop, so it gets updated after running StokesFlow
        self.set_iterative_props(propnames=self.settings['conductance'])

    def setup(self, phase=None, quantity='', conductance='',
              diffusive_conductance='', hydraulic_conductance='', pressure='',
              s_scheme='', **kwargs):
        r"""

        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        if diffusive_conductance:
            self.settings['diffusive_conductance'] = diffusive_conductance
        if hydraulic_conductance:
            self.settings['hydraulic_conductance'] = hydraulic_conductance
        if pressure:
            self.settings['pressure'] = pressure
        if s_scheme:
            self.settings['s_scheme'] = s_scheme
        super().setup(**kwargs)

    def set_outflow_BC(self, pores, mode='merge'):
        r"""
        Adds outflow boundary condition to the selected pores.

        Outflow condition simply means that the gradient of the solved
        quantity does not change, i.e. is 0.

        """
        # Hijack the parse_mode function to verify mode/pores argument
        mode = self._parse_mode(mode, allowed=['merge', 'overwrite', 'remove'],
                                single=True)
        pores = self._parse_indices(pores)
        # Calculating A[i,i] values to ensure the outflow condition
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        throats = network.find_neighbor_throats(pores=pores)
        C12 = network['throat.conns'][throats]
        P12 = phase[self.settings['pressure']][C12]
        gh = phase[self.settings['hydraulic_conductance']][throats]
        Q12 = -gh * np.diff(P12, axis=1).squeeze()
        Qp = np.zeros(self.Np)
        np.add.at(Qp, C12[:, 0], -Q12)
        np.add.at(Qp, C12[:, 1], Q12)
        # Store boundary values
        if ('pore.bc_outflow' not in self.keys()) or (mode == 'overwrite'):
            self['pore.bc_outflow'] = np.nan
        self['pore.bc_outflow'][pores] = Qp[pores]

    def _apply_BCs(self):
        # Apply Dirichlet and rate BCs
        ReactiveTransport._apply_BCs(self)
        if 'pore.bc_outflow' not in self.keys():
            return
        # Apply outflow BC
        diag = self.A.diagonal()
        ind = np.isfinite(self['pore.bc_outflow'])
        diag[ind] += self['pore.bc_outflow'][ind]
        self.A.setdiag(diag)
