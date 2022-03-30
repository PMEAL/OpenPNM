import logging
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import Docorator, SettingsAttr
docstr = Docorator()
logger = logging.getLogger(__name__)

__all__ = ['FickianDiffusion']


@docstr.dedent
class FickianDiffusionSettings():
    r"""

    Parameters
    ----------
    %(ReactiveTransportSettings.parameters)s

    """
    prefix = 'fick'
    quantity = 'pore.concentration'
    conductance = 'throat.diffusive_conductance'


@docstr.dedent
class FickianDiffusion(ReactiveTransport):
    r"""
    A class to simulate binary diffusion with reactions

    Parameters
    ----------
    %(ReactiveTransport.parameters)s

    Notes
    -----
    Fickian diffusion in porous materials occurs in the void space, but
    becuase the diffusion is confined to pores it is impacted by the porosity
    and tortuosity of the network.  Thus the total diffusive flux through the
    network is reduced.  This class can be used to simualte diffusion-reaction
    in domains with specified boundary conditions, or it can be used
    to calculate the effective diffusivity of the network by applying
    controlled boundary conditions on opposing faces, calculating the
    diffusion rate, and inverting Fick's first law as follows:

    .. math::

        D_{eff} = N_{A}*L/(A*\Delta C_{A})

    """

    def __init__(self, settings=None, **kwargs):
        self.settings = SettingsAttr(FickianDiffusionSettings, settings)
        super().__init__(settings=self.settings, **kwargs)
