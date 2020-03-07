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
    ion = None
    quantity = 'pore.concentration'
    conductance = 'throat.ad_dif_mig_conductance'


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
        self.setup(ion=ion)

    def setup(self, phase=None, quantity='', conductance='', ion='', **kwargs):
        r"""
        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
            self.setup(ion=ion)
        if conductance:
            self.settings['conductance'] = conductance
            self.setup(ion=ion)
        if ion:
            if not type(ion) is str:
                ion = ion.name
            self.settings['ion'] = ion
            cond_str = self.settings['conductance']
            cond_str = '.'.join(cond_str.split('.')[:2])
            cond_str += ('.' + ion.name)
            self.settings['conductance'] = cond_str
            quan_str = self.settings['quantity']
            quan_str = '.'.join(quan_str.split('.')[:2])
            quan_str += ('.' + ion.name)
            self.settings['quantity'] = cond_str
        super().setup(**kwargs)
