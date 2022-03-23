import logging
import numpy as np
from openpnm.core import Base
logger = logging.getLogger(__name__)


class GenericMetric(Base):
    r"""
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        project = self.project
        self['pore.all'] = np.ones((project.network.Np, ), dtype=bool)
        self['throat.all'] = np.ones((project.network.Nt, ), dtype=bool)
