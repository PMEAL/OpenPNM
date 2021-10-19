from openpnm.core import Base, LegacyMixin, LabelMixin
from openpnm.utils import logging, Docorator
import numpy as np
logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.get_sections(base='GenericAlgorithm', sections=['Parameters'])
@docstr.dedent
class GenericAlgorithm(Base, LegacyMixin, LabelMixin):
    r"""
    Generic class to define the foundation of Algorithms.

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
        self.settings.setdefault('prefix', 'alg')
        self.settings.update(settings)

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
