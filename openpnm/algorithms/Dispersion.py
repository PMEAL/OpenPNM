"""
===============================================================================
module Dispersion: Class for solving advection diffusion
===============================================================================
"""

import numpy as np
from scipy.sparse.csgraph import laplacian
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class Dispersion(ReactiveTransport):
    r'''
    A subclass of GenericTransport to simulate advection diffusion.

    Examples
    --------
    >>>

    '''
    def __init__(self, settings={}, **kwargs):
        def_set = {'quantity': 'pore.concentration',
                   'diffusive_conductance': 'throat.diffusive_conductance',
                   'hydraulic_conductance': 'throat.hydraulic_conductance',
                   'pressure': 'pore.pressure',
                   'gui': {'setup':        {'quantity': '',
                                            'diffusive_conductance': '',
                                            'hydraulic_conductance': '',
                                            'pressure': ''},
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

    def setup(self, phase=None, quantity='', diffusive_conductance='',
              hydraulic_conductance='', pressure='', **kwargs):
        r"""

        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if diffusive_conductance:
            self.settings['diffusive_conductance'] = diffusive_conductance
        if hydraulic_conductance:
            self.settings['hydraulic_conductance'] = hydraulic_conductance
        if pressure:
            self.settings['pressure'] = pressure
        super().setup(**kwargs)

    def _build_A(self, force=False):
        r"""
        """
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        conns = network['throat.conns']

        P = phase[self.settings['pressure']]
        gh = phase[self.settings['hydraulic_conductance']]
        gd = phase[self.settings['diffusive_conductance']]
        gd = np.tile(gd, 2)

        Qij = -gh*np.diff(P[conns], axis=1).squeeze()
        Qij = np.append(Qij, -Qij)
        Peij = Qij/gd

        Peij[(Peij < 1e-10) & (Peij >= 0)] = 1e-10
        Peij[(Peij > -1e-10) & (Peij <= 0)] = -1e-10
        Qij = Peij*gd

        if force:
            self._pure_A = None
        if self._pure_A is None:
            w = -Qij / (1 - np.exp(Peij))
            A = network.create_adjacency_matrix(weights=w)
            A = laplacian(A)
            self._pure_A = A
        self.A = self._pure_A.copy()
