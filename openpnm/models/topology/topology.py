r"""

.. autofunction:: openpnm.models.topology.coordination_number

"""
from openpnm.utils import logging
logger = logging.getLogger(__name__)


def coordination_number(target):
    r"""
    """
    network = target.network
    N = network.num_neighbors(pores=network.Ps, flatten=False)
    vals = N[network.pores(target.name)]
    return vals
