from openpnm.core import Base, LegacyMixin, LabelMixin
from openpnm.utils import logging, Docorator, SettingsData, SettingsAttr
from traits.api import Str
import numpy as np
logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.get_sections(base='GenericAlgorithmSettings', sections=docstr.all_sections)
@docstr.dedent
class GenericAlgorithmSettings(SettingsData):
    r"""

    Parameters
    ----------
    prefix : str
        The prefix to use when generating a name for the algorithm.

    """
    prefix = Str('alg')


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

    def __init__(self, network=None, project=None, settings={}, **kwargs):
        self.settings._update(GenericAlgorithmSettings(), docs=True)
        self.settings._update(settings)
        super().__init__(project=project, network=network, **kwargs)
        project = self.network.project
        if project:
            self['pore.all'] = np.ones(project.network.Np, dtype=bool)
            self['throat.all'] = np.ones(project.network.Nt, dtype=bool)

    def results(self):
        r"""
        """
        raise NotImplementedError("This method must be subclassed")

    def reset(self):
        r"""
        """
        raise NotImplementedError("This method must be subclassed")
