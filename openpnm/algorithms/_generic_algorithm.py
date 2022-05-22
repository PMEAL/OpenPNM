import logging
import numpy as np
from openpnm.core import Domain
from openpnm.utils import Docorator, SettingsAttr


__all__ = ['GenericAlgorithm']


logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.get_sections(base='GenericAlgorithmSettings', sections=docstr.all_sections)
@docstr.dedent
class GenericAlgorithmSettings:
    r"""

    Parameters
    ----------
    %(BaseSettings.parameters)s

    """


@docstr.get_sections(base='GenericAlgorithm', sections=['Parameters'])
@docstr.dedent
class GenericAlgorithm(Domain):
    r"""
    Generic class to define the foundation of Algorithms

    Parameters
    ----------
    network : OpenPNM Network
        The network object to which this algorithm will apply
    name : str, optional
        Name of the algorithm. If not provided one is generated.

    """

    def __init__(self, network, settings=None, **kwargs):
        self._settings = SettingsAttr(GenericAlgorithmSettings)
        self.settings._update(settings)
        if 'name' not in kwargs.keys():
            kwargs['name'] = 'alg_01'
        super().__init__(network=network, settings=self.settings, **kwargs)
        self['pore.all'] = np.ones(network.Np, dtype=bool)
        self['throat.all'] = np.ones(network.Nt, dtype=bool)
