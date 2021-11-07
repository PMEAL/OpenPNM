from openpnm.algorithms import TransientReactiveTransport, NernstPlanck
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class TransientNernstPlanckSettings:
    phase = None
    quantity = ''
    conductance = ''
    diffusive_conductance = ''


class TransientNernstPlanck(TransientReactiveTransport, NernstPlanck):
    r"""
    A subclass of GenericTransport to perform steady and transient simulations
    of pure diffusion, advection-diffusion and advection-diffusion with
    migration.

    """

    def __init__(self, settings={}, phase=None, ion='', **kwargs):
        super().__init__(**kwargs)
        self.settings._update(TransientNernstPlanckSettings)
        self.settings._update(settings)  # Add user supplied settings
        if phase is not None:
            self.settings['phase'] = phase.name
