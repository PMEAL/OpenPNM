from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


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
        self.setup(ion=ion)

    def setup(self, phase=None, quantity='', conductance='', ion='', **kwargs):
        r"""
        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            if ion:
                if not type(ion) is str:  # Convert ion object to str
                    ion = ion.name
            else:  # Get ion name from settings if not given
                ion = self.settings['ion']
            if not quantity.startswith('pore'):
                quantity = 'pore.' + quantity
            # Remove ion name whether present or not
            quantity = '.'.join(quantity.split('.')[:2])
            quantity += ('.' + ion)  # Re-add ion name
            self.settings['quantity'] = quantity  # Add full value to settings
        if conductance:
            if ion:
                if not type(ion) is str:  # Convert ion object to str
                    ion = ion.name
            else:  # Get ion name from settings if not given
                ion = self.settings['ion']
            if not conductance.startswith('throat'):
                conductance = 'throat.' + conductance
            # Remove ion name whether present or not
            conductance = '.'.join(conductance.split('.')[:2])
            # conductance += ('.' + ion)  # Re-add ion name
            self.settings['conductance'] = conductance
        if ion:
            if not type(ion) is str:
                ion = ion.name
            self.settings['ion'] = ion
            if conductance == '':
                conductance = self.settings['conductance']
            if quantity == '':
                quantity = self.settings['quantity']
            conductance = '.'.join(conductance.split('.')[:2])
            # conductance += ('.' + ion)
            self.settings['conductance'] = conductance
            quantity = '.'.join(quantity.split('.')[:2])
            quantity += ('.' + ion)
            self.settings['quantity'] = quantity
        super().setup(**kwargs)
