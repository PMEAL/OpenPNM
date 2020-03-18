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


class NernstPlanck(ReactiveTransport):
    r"""
    A class to simulate transport of charged species (such as ions) in dilute
    solutions.

    """

    def __init__(self, ion, settings={}, **kwargs):
        super().__init__(**kwargs)
        # self.name = electrolyte  # This interfers with component name
        self.settings._update_settings_and_docs(NernstPlanckSettings())
        self.settings.update(settings)
        # Parse the given ion and append name to quantity and conductance
        if ion:
            if not type(ion) is str:  # Convert ion object to str
                ion = ion.name
            self.settings['ion'] = ion
        quantity = self.settings['quantity']
        if not quantity.endswith(ion):
            quantity = '.'.join(quantity.split('.')[:2])
            quantity += ('.' + ion)  # Re-add ion name
            self.settings['quantity'] = quantity  # Add full value to settings
        conductance = self.settings['conductance']
        if not conductance.endswith(ion):
            conductance = '.'.join(conductance.split('.')[:2])
            conductance += ('.' + ion)  # Re-add ion name
            self.settings['conductance'] = conductance

    def setup(self, phase=None, quantity='', conductance='', ion='', **kwargs):
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
        if ion:
            self.settings['quantity'] = quantity
        super().setup(**kwargs)
