from openpnm.algorithms import AdvectionDiffusion
from openpnm.algorithms import TransientTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class TransientAdvectionDiffusion(AdvectionDiffusion, TransientTransport):
    r"""
    A subclass of GenericTransport to perform steady and transient simulations
    of pure diffusion and advection diffusion problems.

    """

    def __init__(self, t_initial=0, t_final=8000, t_step=2, t_output=200,
                 tolerance=1e-4,
                 t_scheme='implicit', s_scheme='powerlaw',
                 settings={}, **kwargs):
        super().__init__(**kwargs)
        # Set some default settings
        self.settings.update({
                'quantity': 'pore.mole_fraction',
                'diffusive_conductance': 'throat.diffusive_conductance',
                'hydraulic_conductance': 'throat.hydraulic_conductance',
                'pressure': 'pore.pressure',
                'molar_density': 'pore.molar_density',
                't_initial': t_initial,
                't_final': t_final,
                't_step': t_step,
                't_output': t_output,
                'tolerance': tolerance,
                't_scheme': t_scheme,
                's_scheme': s_scheme})
        # Apply any received settings to overwrite defaults
        self.settings.update(settings)
        # Save A matrix of the steady sys of eqs
        self.A = self.build_A()
