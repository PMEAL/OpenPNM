import numpy as np
from openpnm.core import Base, LegacyMixin, LabelMixin
from openpnm.utils import logging, Docorator, SettingsAttr
logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.get_sections(base='GenericAlgorithmSettings', sections=docstr.all_sections)
@docstr.dedent
class GenericAlgorithmSettings:
    r"""

    Parameters
    ----------
    prefix : str
        The prefix to use when generating a name for the algorithm.

    """
    prefix = 'alg'
    name = ''


@docstr.get_sections(base='GenericAlgorithm', sections=['Parameters'])
@docstr.dedent
class GenericAlgorithm(Base, LegacyMixin, LabelMixin):
    r"""
    Generic class to define the foundation of Algorithms

    Parameters
    ----------
    network : GenericNetwork
        The network object to which this algorithm will apply.
    name : str, optional
        Name of the algorithm
    project : Project, optional
        Either a Network or a Project must be supplied

    """

    def __init__(self, settings={}, **kwargs):
        self.settings = SettingsAttr(GenericAlgorithmSettings, settings)
        super().__init__(settings=self.settings, **kwargs)

        if self.project:
            self['pore.all'] = np.ones(self.project.network.Np, dtype=bool)
            self['throat.all'] = np.ones(self.project.network.Nt, dtype=bool)

    def results(self):
        r"""
        """
        raise NotImplementedError("This method must be subclassed")

    def reset(self):
        r"""
        """
        raise NotImplementedError("This method must be subclassed")
