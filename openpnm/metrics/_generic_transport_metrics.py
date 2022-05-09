import logging
import numpy as np
from openpnm.metrics import GenericMetric
from openpnm.topotools import get_domain_area, get_domain_length
logger = logging.getLogger(__name__)


class GenericTransportMetrics(GenericMetric):
    r"""
    """
    def _calc_eff_prop(self, inlets=None, outlets=None,
                       domain_area=None, domain_length=None, rates=None,
                       prop_diff=None):
        r"""
        Calculates the effective transport through the network.
        Parameters
        ----------
        inlets : array_like
            The pores where the inlet boundary conditions were applied. If
            not given an attempt is made to infer them from the algorithm.
        outlets : array_like
            The pores where the outlet boundary conditions were applied.
            If not given an attempt is made to infer them from the
            algorithm.
        domain_area : float
            The area of the inlet and/or outlet face (which shold match)
        domain_length : float
            The length of the domain between the inlet and outlet faces
        Returns
        -------
        The effective transport property through the network
        """
        network = self.project.network
        if rates is None:
            raise Exception('You need to run the algorithm first.')
        Dx = prop_diff
        if domain_area is None:
            domain_area = get_domain_area(network, inlets=inlets, outlets=outlets)
        if domain_length is None:
            domain_length = get_domain_length(network, inlets=inlets, outlets=outlets)
        D = np.sum(rates) * domain_length / domain_area / Dx
        return D
