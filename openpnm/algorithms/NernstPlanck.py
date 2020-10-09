import numpy as np
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sectionsf('NernstPlanckSettings', sections=['Parameters'])
@docstr.dedent
class NernstPlanckSettings(GenericSettings):
    r"""
    Parameters
    ----------
    %(ReactiveTransportSettings.parameters)s
    quantity : string (default = 'pore.concentration')
        The quantity to solve for.  Note that this will have the 'ion' name
        appended to the end (i.e. ``'pore.concentration.Na'``)
    conductance : string (default is 'throat.ad_dif_mig_conductance')
        The conductance of the ion.

    Other Parameters
    ----------------
    s_scheme : string (default = 'exponential')
        ##

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
    diffusive_conductance = 'throat.diffusive_conductance'


class NernstPlanck(ReactiveTransport):
    r"""
    A class to simulate transport of charged species (such as ions) in dilute
    solutions.

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
        gd = phase[self.settings['diffusive_conductance']]
        gd = gd[throats]
        gt = phase[self.settings['conductance']][:, 0]
        gt = gt[throats]
        # remove the diffusive contribution (only keep adv and mig)
        Q12 = gd-gt
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
