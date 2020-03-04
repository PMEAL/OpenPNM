from openpnm.algorithms import TransientReactiveTransport, FickianDiffusion
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class TransientFickianDiffusion(TransientReactiveTransport, FickianDiffusion):
    r"""
    A class to simulate transient diffusion with reactions

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        super().__init__(**kwargs)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    def setup(self, phase=None, quantity='', conductance='',
              t_initial=None, t_final=None, t_step=None, t_output=None,
              t_tolerance=None, t_scheme='', **kwargs):
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        if t_initial:
            self.settings['t_initial'] = t_initial
        if t_final:
            self.settings['t_final'] = t_final
        if t_step:
            self.settings['t_step'] = t_step
        if t_output is not None:
            self.settings['t_output'] = t_output
        if t_tolerance:
            self.settings['t_tolerance'] = t_tolerance
        if t_scheme:
            self.settings['t_scheme'] = t_scheme
        self.settings.update(kwargs)
