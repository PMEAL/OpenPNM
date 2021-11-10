from openpnm.algorithms import TransientReactiveTransport, NernstPlanck
from openpnm.utils import logging, SettingsAttr
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
        self.settings = SettingsAttr(TransientNernstPlanckSettings, settings)
        super().__init__(settings=self.settings, **kwargs)
        if phase is not None:
            self.settings['phase'] = phase.name
