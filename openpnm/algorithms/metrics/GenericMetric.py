from openpnm.core import Base
from openpnm.utils import logging
import numpy as np
logger = logging.getLogger(__name__)


class GenericMetric(Base):
    r"""


    """

    def __init__(self, network=None, project=None, settings={}, **kwargs):
        self.settings.setdefault('prefix', 'alg')
        self.settings.update(settings)
        if network is not None:
            project = network.project
        super().__init__(project=project, **kwargs)

        # Deal with network or project arguments
        if network is not None:
            if project is not None:
                assert network is project.network
            else:
                project = network.project
        if project:
            self['pore.all'] = np.ones((project.network.Np, ), dtype=bool)
            self['throat.all'] = np.ones((project.network.Nt, ), dtype=bool)
