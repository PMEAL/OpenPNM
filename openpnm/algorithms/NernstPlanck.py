from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class NernstPlanck(ReactiveTransport):
    r"""
    A class to simulate transport of charged species (such as ions) in dilute
    solutions.

    """
    def __init__(self, settings={}, phase=None, ion='', **kwargs):
        def_set = {'phase': None,
                   'quantity': 'pore.concentration.'+ion,
                   'conductance': 'throat.ad_dif_mig_conductance.'+ion,
                   'ion': ion}
        super().__init__(**kwargs)
        # self.name = electrolyte  # This interfers with component name
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    def setup(self, phase=None, quantity='', conductance='', ion='', **kwargs):
        r"""
        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        if ion:
            self.settings['ion'] = ion
        super().setup(**kwargs)
