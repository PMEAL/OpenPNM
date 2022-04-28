import logging
import numpy as np
from openpnm.core import Base, LegacyMixin, Base2, Domain
from openpnm.utils import Docorator, SettingsAttr
logger = logging.getLogger(__name__)
docstr = Docorator()

__all__ = ['GenericAlgorithm']


@docstr.get_sections(base='GenericAlgorithmSettings', sections=docstr.all_sections)
@docstr.dedent
class GenericAlgorithmSettings:
    r"""

    Parameters
    ----------
    %(BaseSettings.parameters)s

    """
    prefix = 'alg'
    name = ''


@docstr.get_sections(base='GenericAlgorithm', sections=['Parameters'])
@docstr.dedent
class GenericAlgorithm(Domain, LegacyMixin):
    r"""
    Generic class to define the foundation of Algorithms

    Parameters
    ----------
    network : GenericNetwork
        The network object to which this algorithm will apply.
    name : str, optional
        Name of the algorithm

    """

    def __init__(self, network, settings=None, **kwargs):
        self.settings = SettingsAttr(GenericAlgorithmSettings, settings)
        super().__init__(network=network, settings=self.settings, **kwargs)
        self['pore.all'] = np.ones(network.Np, dtype=bool)
        self['throat.all'] = np.ones(network.Nt, dtype=bool)

    def results(self):
        r"""
        """
        raise NotImplementedError("This method must be subclassed")

    def reset(self):
        r"""
        """
        raise NotImplementedError("This method must be subclassed")
