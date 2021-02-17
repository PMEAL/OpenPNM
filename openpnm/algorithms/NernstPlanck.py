import numpy as np
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sections(base='NernstPlanckSettings', sections=['Parameters'])
@docstr.dedent
class NernstPlanckSettings(GenericSettings):
    r"""
    Parameters
    ----------
    %(ReactiveTransportSettings.parameters)s
    quantity : str, optional
        The quantity to solve for. The default value is
        'pore.concentration'. Note that this will have the 'ion' name
        appended to the end (i.e. ``pore.concentration.Na``)
    conductance : str, optional
        The overall mass transfer conductance of the ion. The default
        value is 'throat.ad_dif_mig_conductance'.
    diffusive_conductance : str
        The name of the diffusive conductance values. The default value is
        'throat.diffusive_conductance'.
    hydraulic_conductance : str, optional
        The name of the hydraulic conductance values. The default value is
        'throat.hydraulic_conductance'.
    pressure : str, optional
        The name of the pressure values calculated by the ``StokesFlow``
        algorithm. The default value is 'pore.pressure'.

    Other Parameters
    ----------------
    s_scheme : str
        The discretization scheme for the advective terms. The default
        value is 'powerlaw'.

    ----

    **The following parameters pertain to the ReactiveTransport class**

    %(ReactiveTransportSettings.other_parameters)s

    ----

    **The following parameters pertain to the GenericTransport class**

    %(GenericTransportSettings.other_parameters)s

    """
    ion = ''
    quantity = 'pore.concentration'
    conductance = 'throat.ad_dif_mig_conductance'
    hydraulic_conductance = 'throat.hydraulic_conductance'
    diffusive_conductance = 'throat.diffusive_conductance'
    pressure = 'pore.pressure'


class NernstPlanck(ReactiveTransport):
    r"""
    A subclass of ``ReactiveTransport`` to simulate transport of charged
    species (such as ions) in dilute solutions.

    """

    def __init__(self, ion=None, settings={}, **kwargs):
        super().__init__(**kwargs)
        # self.name = electrolyte  # This interfers with component name
        self.settings._update_settings_and_docs(NernstPlanckSettings())
        self.settings.update(settings)
        # Parse the given ion and append name to quantity and conductance
        if ion:
            if not isinstance(ion, str):  # Convert ion object to str
                ion = ion.name
            self.settings['ion'] = ion
        quantity = self.settings['quantity']
        if (ion is not None and not quantity.endswith(ion)):
            quantity = '.'.join(quantity.split('.')[:2])
            quantity += ('.' + ion)  # Re-add ion name
            self.settings['quantity'] = quantity  # Add full value to settings
        conductance = self.settings['conductance']
        if (ion is not None and not conductance.endswith(ion)):
            conductance = '.'.join(conductance.split('.')[:2])
            conductance += ('.' + ion)  # Re-add ion name
            self.settings['conductance'] = conductance
        diffusive_conductance = self.settings['diffusive_conductance']
        if (ion is not None and not diffusive_conductance.endswith(ion)):
            diffusive_conductance = '.'.join(
                diffusive_conductance.split('.')[:2])
            diffusive_conductance += ('.' + ion)  # Re-add ion name
            self.settings['diffusive_conductance'] = diffusive_conductance

    def setup(self, phase=None, quantity='', conductance='',
              diffusive_conductance='', ion='', **kwargs):
        r"""

        Parameters
        ----------
        %(NernstPlanckSettings.parameters)s

        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity  # Add full value to settings
        if conductance:
            self.settings['conductance'] = conductance
        if diffusive_conductance:
            self.settings['diffusive_conductance'] = diffusive_conductance
        if ion:
            self.settings['quantity'] = quantity
        super().setup(**kwargs)

    def set_outflow_BC(self, pores, mode='merge'):
        r"""
        Adds outflow boundary condition to the selected pores.

        Notes
        -----
        Outflow condition simply means that the gradient of the solved
        quantity does not change, i.e. is 0.

        """
        # Hijack the parse_mode function to verify mode/pores argument
        allowed_modes = ['merge', 'overwrite', 'remove']
        mode = self._parse_mode(mode, allowed=allowed_modes, single=True)
        pores = self._parse_indices(pores)

        # Calculating A[i,i] values to ensure the outflow condition
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        throats = network.find_neighbor_throats(pores=pores)
        C12 = network.conns[throats]
        P12 = phase[self.settings['pressure']][C12]
        # Note: only keep the advection term when flow leaves the outflow pores
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
