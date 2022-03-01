import numpy as np
from openpnm.core import Base
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class GenericMetric(Base):
    """Brief description of 'GenericMetric'"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        project = self.project
        self['pore.all'] = np.ones((project.network.Np, ), dtype=bool)
        self['throat.all'] = np.ones((project.network.Nt, ), dtype=bool)
