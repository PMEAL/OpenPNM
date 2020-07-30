from openpnm.core import Base
from openpnm.utils import logging
import numpy as np
logger = logging.getLogger(__name__)


class GenericMetric(Base):
    r"""

    """

    def __init__(self, network=None, project=None, settings={}, **kwargs):
        self.settings['prefix'] = 'metric'
        super().__init__(project=project, **kwargs)
        project = self.project
        self['pore.all'] = np.ones((project.network.Np, ), dtype=bool)
        self['throat.all'] = np.ones((project.network.Nt, ), dtype=bool)
