import numpy as np
from scipy.sparse.csgraph import laplacian
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class AdvectionDiffusion(ReactiveTransport):
    r"""
    A subclass of GenericTransport to simulate advection diffusion

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'quantity': 'pore.concentration',
                   'diffusive_conductance': 'throat.diffusive_conductance',
                   'hydraulic_conductance': 'throat.hydraulic_conductance',
                   'pressure': 'pore.pressure',
                   's_scheme': 'powerlaw',
                   'gui': {'setup':        {'phase': None,
                                            'quantity': '',
                                            'diffusive_conductance': '',
                                            'hydraulic_conductance': '',
                                            'pressure': '',
                                            's_scheme': ''},
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

    def setup(self, phase=None, quantity='', diffusive_conductance='',
              hydraulic_conductance='', pressure='', s_scheme='', **kwargs):
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
        if s_scheme:
            self.settings['s_scheme'] = s_scheme
        super().setup(**kwargs)

    def _build_A(self, force=False):
        s_dis = self.settings['s_scheme']
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
            if s_dis == 'upwind':
                w = gd + np.maximum(0, -Qij)
                A = network.create_adjacency_matrix(weights=w)
            elif s_dis == 'hybrid':
                w = np.maximum(0, np.maximum(-Qij, gd-Qij/2))
                A = network.create_adjacency_matrix(weights=w)
            elif s_dis == 'powerlaw':
                w = gd * np.maximum(0, (1 - 0.1*np.abs(Peij))**5) + \
                    np.maximum(0, -Qij)
                A = network.create_adjacency_matrix(weights=w)
            elif s_dis == 'exponential':
                w = -Qij / (1 - np.exp(Peij))
                A = network.create_adjacency_matrix(weights=w)
            else:
                raise Exception('Unrecognized discretization scheme: ' + s_dis)
            A = laplacian(A)
            self._pure_A = A
        self.A = self._pure_A.copy()

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

    def rate(self, pores=[], throats=[], mode='group'):
        r"""
        Calculates the net rate of material moving into a given set of pores or
        throats

        Parameters
        ----------
        pores : array_like
            The pores for which the rate should be calculated

        throats : array_like
            The throats through which the rate should be calculated

        mode : string, optional
            Controls how to return the rate.  Options are:

            *'group'*: (default) Returns the cumulative rate of material
            moving into the given set of pores

            *'single'* : Calculates the rate for each pore individually

        Returns
        -------
        If ``pores`` are specified, then the returned values indicate the
        net rate of material exiting the pore or pores.  Thus a positive
        rate indicates material is leaving the pores, and negative values
        mean material is entering.

        If ``throats`` are specified the rate is calculated in the direction of
        the gradient, thus is always positive.

        If ``mode`` is 'single' then the cumulative rate through the given
        pores (or throats) are returned as a vector, if ``mode`` is 'group'
        then the individual rates are summed and returned as a scalar.

        """
        pores = self._parse_indices(pores)
        throats = self._parse_indices(throats)

        network = self.project.network
        C12 = network["throat.conns"]
        phase = self.project.phases()[self.settings['phase']]

        P = phase[self.settings['pressure']]
        X = self[self.settings['quantity']]
        gh = phase[self.settings['hydraulic_conductance']]
        gd = phase[self.settings['diffusive_conductance']]

        Q12 = -gh * np.diff(P[C12], axis=1).squeeze()
        Pe12 = Q12 / gd
        Pe12 = np.abs(Pe12).clip(min=1e-10) * np.sign(Pe12)
        Q12 = Pe12 * gd

        X12 = X[C12]
        dX = -np.diff(X12, axis=1).squeeze()

        s_dis = self.settings['s_scheme']

        if s_dis == 'upwind':
            w = gd + np.maximum(0, -Q12)
            Qt = -Q12 * X12[:, 0] - dX * w
        elif s_dis == 'hybrid':
            w = np.maximum(0, np.maximum(-Q12, gd-Q12/2))
            Qt = -Q12 * X12[:, 0] - dX * w
        elif s_dis == 'powerlaw':
            w = gd * np.maximum(0, (1 - 0.1*np.abs(Pe12))**5) + \
                np.maximum(0, -Q12)
            Qt = -Q12 * X12[:, 0] - dX * w
        elif s_dis == 'exponential':
            w = -Q12 / (1 - np.exp(Pe12))
            Qt = -Q12 * X12[:, 0] - dX * w

        if len(throats) and len(pores):
            raise Exception('Must specify either pores or throats, not both')

        if len(throats):
            R = np.absolute(Qt[throats])

        if len(pores):
            Qp = np.zeros((self.Np,))
            np.add.at(Qp, C12[:, 0], -Qt)
            np.add.at(Qp, C12[:, 1], Qt)
            R = Qp[pores]

        if mode == 'group':
            R = np.sum(R)

        return np.array(R, ndmin=1)
