from openpnm.algorithms import StokesFlow
from openpnm.algorithms import TransientTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class TransientStokesFlow(StokesFlow, TransientTransport):
    r"""
    A subclass of GenericTransport to simulate Stokes flow.

    """

    def __init__(self, t_initial=0, t_final=1e+06, t_step=2, t_output=50,
                 tolerance=1e-5, t_scheme='cranknicolson',
                 settings={}, **kwargs):
        super().__init__(**kwargs)
        # Set some default settings
        self.settings.update({'quantity': 'pore.pressure',
                              'conductance': 'throat.hydraulic_conductance',
                              't_initial': t_initial,
                              't_final': t_final,
                              't_step': t_step,
                              't_output': t_output,
                              'tolerance': tolerance,
                              't_scheme': t_scheme})
        # Apply any received settings to overwrite defaults
        self.settings.update(settings)
        # Save A matrix of the steady sys of eqs
        self.A = self.build_A()
        # Define _coef as unity
        self._coef = 1
