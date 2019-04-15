from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class NernstPlanck(ReactiveTransport):
    r"""
    """
    def __init__(self, settings={}, phase=None, electrolyte='', **kwargs):
        def_set = {'phase': None,
                   'quantity': 'pore.concentration.'+electrolyte,
                   'conductance': 'throat.ad_dif_mig_conductance.'+electrolyte}
        super().__init__(**kwargs)
        self.name = electrolyte
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    def setup(self, phase=None, quantity='', conductance='', **kwargs):
        r"""
        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        super().setup(**kwargs)
